"""
    Electrum.PlanewaveIndex{D}

A special indexing type used to index the components of wavefunctions in a planewave basis.

In many computational chemistry packages, the standard indexing of wavefunction components occurs in
the following canonical order: spins, k-points, bands, then the h, k, and l indices of the G-vectors
associated with the coefficients. However, in many cases, users will want to select a spin, k-point,
or band before selecting a G-vector index.

To keep the syntax intuitive while maintaining performance (mostly by ensuring that all of the
components of a wavefunction are stored compactly) this index type ensures that the iteration occurs
in a natural order for users and in an efficient order for Julia.
    
This also ensures that G-vectors with negative indices are handled correctly: while the canonical
G-vectors are those within some defined ranges of indices, out of bounds indices are reinterpreted
automatically using modulo arithmetic.
"""
struct PlanewaveIndex{D}
    spin::Int
    kpoint::Int
    band::Int
    g::CartesianIndex{D}
    PlanewaveIndex(spin, kpt, band, g::CartesianIndex{D}) where D = new{D}(spin, kpt, band, g)
end

PlanewaveIndex(spin, kpt, band, g...) = PlanewaveIndex(spin, kpt, band, CartesianIndex(g))
PlanewaveIndex(spin, kpt, band, g::NTuple) = PlanewaveIndex(spin, kpt, band, CartesianIndex(g))
PlanewaveIndex(spin, kpt, band, g::StaticVector) = PlanewaveIndex(spin, kpt, band, Tuple(g))
PlanewaveIndex(i::CartesianIndex) = PlanewaveIndex(Tuple(i)...)

function PlanewaveIndex{D}(spin, kpt, band, gvec::AbstractVector) where D
    return PlanewaveIndex(spin, kpt, band, ntuple(i -> gvec[i], Val{D}()))
end

Base.show(io::IO, i::PlanewaveIndex) = print(io, PlanewaveIndex, (i.spin, i.kpoint, i.band, i.g))

Tuple(i::PlanewaveIndex) = (i.spin, i.kpoint, i.band, Tuple(i.g)...)
CartesianIndex(i::PlanewaveIndex) = CartesianIndex(Tuple(i))

Base.convert(T::Type{<:Tuple}, i::PlanewaveIndex) = Tuple(i)::T
Base.convert(T::Type{<:CartesianIndex}, i::PlanewaveIndex) = CartesianIndex(i)::T
Base.convert(T::Type{<:PlanewaveIndex}, i::CartesianIndex) = PlanewaveIndex(i)::T

#=
"""
    Electrum.PlanewaveIndices{D} <: AbstractArray{PlanewaveIndex{D},D}

An array of valid indices of a `PlanewaveWavefunction{D,T}` within the G-vector bounds determined by
the energy cutoff of the calculation that generated the data.
"""
struct PlanewaveIndices{D} <: AbstractArray{PlanewaveIndex{D},D}
    spins::Base.OneTo{Int}
    kpoints::Base.OneTo{Int}
    bands::Base.OneTo{Int}
    grange::NTuple{D,UnitRange{Int}}
    function PlanewaveIndices(
        spins::AbstractUnitRange,
        kpoints::AbstractUnitRange,
        bands::AbstractUnitRange,
        grange::NTuple{D,<:AbstractUnitRange}
    ) where D
        return new{D}(spins, kpoints, bands, grange)
    end
end

function PlanewaveIndices(
    spins::AbstractUnitRange,
    kpts::AbstractUnitRange,
    bands::AbstractUnitRange,
    gs::AbstractUnitRange...
)
    return PlanewaveIndices(spins, kpts, bands, gs)
end

function PlanewaveIndices(spins::Integer, kpts::Integer, bands::Integer, gs::AbstractUnitRange...)
    return PlanewaveIndices(Base.OneTo(spins), Base.OneTo(kpts), Base.OneTo(bands), gs)
end

Base.axes(p::PlanewaveIndices) = (p.spins, p.kpoints, p.bands, p.grange...)
Base.size(p::PlanewaveIndices) = length.(axes(p))

function Base.show(io::IO, p::PlanewaveIndices)
    print(
        io, PlanewaveIndices,
        (last(p.spins), last(p.kpoints), last(p.bands), p.grange...)
    )
end

Base.show(io::IO, ::MIME"text/plain", p::PlanewaveIndices) = show(io, p)

# Convert to a vector within the allowed range of G-vectors
# In other words, favor a negative index in that range over a larger positive index
function Base.getindex(p::PlanewaveIndices{D}, i::PlanewaveIndex{D}) where D
    i.spin in p.spins || throw(BoundsError(p, i))
    i.kpoint in p.kpoints || throw(BoundsError(p, i))
    i.band in p.bands || throw(BoundsError(p, i))
    sz = length.(p.grange)
    c = CartesianIndex(mod.(Tuple(i.g) .+ div.(sz, 2), sz) .- div.(sz, 2))
    return PlanewaveIndex(i.spin, i.kpoint, i.band, c)
end

Base.getindex(p::PlanewaveIndices, ::Integer, ::Integer, ::Integer) = CartesianIndices(p.grange)

function Base.getindex(
    p::PlanewaveIndices{D},
    spin::Integer,
    kpt::Integer,
    band::Integer,
    i::Vararg{<:Integer,D}    
) where D
    return p[PlanewaveIndex(spin, kpt, band, i...)]
end

Base.LinearIndices(p::PlanewaveIndices) = LinearIndices((p.grange..., p.bands, p.kpoints, p.spins))
=#

"""
    PlanewaveWavefunction{D,T} <: AbstractArray{T,D}

Stores the components of a wavefunction constructed from a planewave basis. Usually, the coefficient
data type `T` will be a `ComplexF32`, as in DFT calculations, double precision convergence of the
density will correspond to single-precision converegnce of the wavefunction.

Internally, coefficients are stored in an `Array{4,T}`. Indexing is then manually implemented, with
a `D`-dimensional `CartesianIndex` used for accessing each coefficient associated with a G-vector.
`PlanewaveWavefunction` instances are mutable, with `getindex()` and `setindex!()` defined for them,
but they are not resizable, and the backing array should not be resized.
"""
struct PlanewaveWavefunction{D,T} <: AbstractArray{T,D}
    basis::ReciprocalBasis{D,Float64}
    spins::Vector{SVector{D,Float64}}
    kpoints::KPointMesh{D}
    energies::Array{Float64,3}
    occupancies::Array{Float64,3}
    grange::NTuple{D,UnitRange{Int}}
    data::Array{T,4}
    function PlanewaveWavefunction(
        basis::LatticeBasis,
        spins::AbstractVector{<:StaticVector{D,<:Real}},
        kpoints::AbstractVector{KPoint{D}},
        energies::AbstractArray{<:Real,3},      # Not specific to 3D, this is from kpt/band/spin
        occupancies::AbstractArray{<:Real,3},
        grange::NTuple{D,<:AbstractUnitRange{<:Integer}},
        data::Array{T,4}
    ) where {D,T}
        @assert length(spins) === size(data, 4) "Mismatch in the number of spins"
        @assert length(kpoints) === size(data, 3) "Mismatch in the number of k-points"
        @assert size(energies) === size(data)[2:4] "Mismatch in the size of energy data"
        @assert size(occupancies) === size(data)[2:4] "Mismatch in the size of occupancy data"
        @assert prod(length.(grange)) === size(data, 1) "G-vector limits do not match data size"
        return new{D,T}(basis, spins, kpoints, energies, occupancies, grange, data)
    end
end

Base.has_offset_axes(::PlanewaveWavefunction) = true

"""
    PlanewaveWavefunction{D,T}(
        basis::LatticeBasis,
        nspin::Integer,
        nkpt::AbstractVector,
        nband::Integer,
        grange::AbstractUnitRange{<:Integer}...
    )

Constructs an empty `PlanewaveWavefunction` with `nspin` spins, a list of k-points `kptmesh`, 
`nband` bands, and G-vectors in the ranges given by `grange`.
"""
function PlanewaveWavefunction{D,T}(
    basis::LatticeBasis,
    nspin::Integer,
    kptmesh::AbstractVector{<:AbstractVector},
    nband::Integer,
    grange::Vararg{AbstractUnitRange{<:Integer},D}
) where {D,T}
    nkpt = length(kptmesh)
    return PlanewaveWavefunction(
        basis,
        zeros(SVector{D,Float64}, nspin),
        kptmesh,
        zeros(Float64, nband, nkpt, nspin),
        zeros(Float64, nband, nkpt, nspin),
        grange,
        zeros(T, prod(length.(grange)), nband, nkpt, nspin)
    )
end

"""
    PlanewaveWavefunction{D,T}(
        basis::LatticeBasis,
        nspin::Integer,
        nkpt::Integer,
        nband::Integer,
        grange::AbstractUnitRange{<:Integer}...
    )

Constructs an empty `PlanewaveWavefunction` with `nspin` spins, `nkpt` k-points, `nband` bands, and
G-vectors in the ranges given by `grange`.
"""
function PlanewaveWavefunction{D,T}(
    basis::LatticeBasis,
    nspin::Integer,
    nkpt::Integer,
    nband::Integer,
    grange::Vararg{AbstractUnitRange{<:Integer},D}
) where {D,T}
    return PlanewaveWavefunction(
        basis,
        zeros(SVector{D,Float64}, nspin),
        KPointMesh(zeros(KPoint{D}, nkpt)),
        zeros(Float64, nband, nkpt, nspin),
        zeros(Float64, nband, nkpt, nspin),
        grange,
        zeros(T, prod(length.(grange)), nband, nkpt, nspin)
    )
end

Base.size(wf::PlanewaveWavefunction) = (reverse(size(wf.energies))..., length.(wf.grange)...)
Base.axes(wf::PlanewaveWavefunction) = (reverse(axes(wf.energies))..., wf.grange...)

# PlanewaveIndices(wf::PlanewaveWavefunction) = PlanewaveIndices(axes(wf)...)

function Base.LinearIndices(wf::PlanewaveWavefunction)
    return LinearIndices((wf.grange..., wf.bands, wf.kpoints, wf.spins))
end

"""
    nspin(wf::PlanewaveWavefunction) -> Int

Returns the number of spins associated with a `PlanewaveWavefunction`.
"""
nspin(wf::PlanewaveWavefunction) = length(wf.spins)

"""
    nband(wf::PlanewaveWavefunction) -> Int

Returns the number of bands (occupied and unoccupied) associated with a `PlanewaveWavefunction`.
"""
nband(wf::PlanewaveWavefunction) = size(wf.energies, 1)

KPointMesh(wf::PlanewaveWavefunction) = wf.kpoints

# Override some of the more generic AbstractDataGrid methods
Base.getindex(wf::PlanewaveWavefunction, i...) = throw(MethodError(getindex, (wf, i...)))
Base.setindex!(wf::PlanewaveWavefunction, i...) = throw(MethodError(setindex!, (wf, i...)))

function Base.getindex(wf::PlanewaveWavefunction{D}, i::PlanewaveIndex{D}) where D
    all(in.((i.spin, i.kpoint, i.band), axes(wf)[1:3])) || throw(BoundsError(wf, Tuple(i)))
    l = LinearIndices(wf.grange)[CartesianIndex(mod.(Tuple(i.g), length.(wf.grange)) .+ 1)]
    @inbounds return wf.data[l, i.band, i.kpoint, i.spin]
end

Base.getindex(wf::PlanewaveWavefunction, i::CartesianIndex) = wf[PlanewaveIndex(i)]

function Base.getindex(wf::PlanewaveWavefunction, spin, kpt, band)
    # Convert integers to ranges so that arrays are sliced into smaller arrays
    s = spin isa Integer ? (spin:spin) : spin
    k = kpt isa Integer ? (kpt:kpt) : kpt
    b = band isa Integer ? (band:band) : band
    return PlanewaveWavefunction(
        wf.basis,
        wf.spins[s],
        wf.kpoints[k],
        wf.energies[b, k, s],
        wf.occupancies[b, k, s],
        wf.grange,
        wf.data[:, b, k, s]
    )
end

function Base.getindex(wf::PlanewaveWavefunction, spin::Integer, kpt::Integer, band::Integer)
    return ReciprocalDataGrid(
        reshape(wf.data[:, band, kpt, spin], length.(wf.grange)),
        basis(wf),
        wf.kpoints[kpt]
    )
end

function Base.getindex(
    wf::PlanewaveWavefunction{D},
    spin::Integer,
    kpt::Integer,
    band::Integer,
    g::Vararg{Integer,D}
) where D
    return wf[PlanewaveIndex(spin, kpt, band, g)]
end

function Base.setindex!(wf::PlanewaveWavefunction{D}, x, i::PlanewaveIndex{D}) where D
    l = LinearIndices(wf.grange)[CartesianIndex(mod.(Tuple(i.g), length.(wf.grange)) .+ 1)]
    wf.data[l, i.band, i.kpoint, i.spin] = x
end

Base.setindex!(w::PlanewaveWavefunction, x, i::CartesianIndex) = setindex!(w, x, PlanewaveIndex(i))

function Base.setindex!(
    wf::PlanewaveWavefunction{D},
    x,
    spin::Integer,
    kpt::Integer,
    band::Integer,
    g::Vararg{Integer,D}
) where D
    wf[PlanewaveIndex(spin, kpt, band, g)] = x
end

function Base.show(io::IO, wf::T) where T<:PlanewaveWavefunction
    print(io, PlanewaveWavefunction, Tuple(getfield(wf, s) for s in fieldnames(T)))
end

FFTBins(wf::PlanewaveWavefunction) = FFTBins(length.(wf.grange)...)

"""
    fermi(wf::PlanewaveWavefunction) -> Float64

Estimates the Fermi energy associated with a reciprocal space wavefunction using the energy and
occupancy data in the `PlanewaveWavefunction`.
"""
function fermi(wf::PlanewaveWavefunction)
    # Generate a matrix of energies, occupancies, and indices
    eo = collect(zip(wf.energies, wf.occupancies))
    # Get the maximum occupancy
    maxocc = round(Int, maximum(x -> x[2], eo))
    @assert maxocc in 1:2 "The calculated maximum occupancy was $maxocc."
    # Convert it to a vector sorted by energies
    eo_sorted = sort!(vec(eo), by=(x -> x[1]))
    # Find the index of the last occupancy that's greater than half of maxocc
    ind = findlast(x -> x[2] > maxocc/2, eo_sorted)
    # Perform a linear interpolation between O[ind] and O[ind+1] 
    approx = (maxocc/2 - eo_sorted[ind][2]) / (eo_sorted[ind+1][2] - eo_sorted[ind][2])
    return eo_sorted[ind][1] * approx + eo_sorted[ind+1][1] * (1 - approx)
end

energies(wf::PlanewaveWavefunction) = wf.energies
occupancies(wf::PlanewaveWavefunction) = wf.occupancies
EnergiesOccupancies(wf::PlanewaveWavefunction) = EnergyOccupancy.(wf.energies, wf.occupancies)

function EnergiesOccupancies{T}(wf::PlanewaveWavefunction) where T
    return EnergyOccupancy{T}.(wf.energies, wf.occpuancies)
end

"""
    min_energy(wf::PlanewaveWavefunction) -> Float64

Returns the minimum energy value (in Hartrees) in the energy entries for a `PlanewaveWavefunction`.
"""
min_energy(wf::PlanewaveWavefunction) = minimum(wf.energies)

"""
    max_energy(wf::PlanewaveWavefunction) -> Float64

Returns the maximum energy value (in Hartrees) in the energy entries for a `PlanewaveWavefunction`.

Note that this maximum energy will likely correspond to an unoccupied state, and should not be taken
to be the Fermi energy. For this value, see `fermi(wf)`.
"""
max_energy(wf::PlanewaveWavefunction) = maximum(wf.energies)

"""
    max_occupancy(wf::PlanewaveWavefunction) -> Int

Returns the maximum occupancy value associated with a `PlanewaveWavefunction`. This value should be
2 if `nspin(wf) === 1`, and 1 otherwise.
"""
max_occupancy(wf::PlanewaveWavefunction) = round(Int, maximum(wf.occupancies))
