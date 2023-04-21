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

function PlanewaveIndex{D}(spin, kpt, band, gvec::AbstractVector) where D
    return PlanewaveIndex(spin, kpt, band, ntuple(i -> gvec[i], Val{D}()))
end

Base.show(io::IO, i::PlanewaveIndex) = print(io, PlanewaveIndex, (i.spin, i.kpoint, i.band, i.g))

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

Base.getindex(p::PlanewaveIndices, spin, kpt, band) = CartesianIndices(p.grange)
Base.LinearIndices(p::PlanewaveIndices) = LinearIndices((p.grange..., p.bands, p.kpoints, p.spins))

function Base.iterate(p::PlanewaveIndices, i = 1)
    (i in 1:length(p)) || return nothing
    # Get the spin, k-point, and band
    band = (i % prod(length.(p.grange))) + 1
end

#---The meat and potatoes--------------------------------------------------------------------------#
"""
    PlanewaveWavefunction{D,T} <: AbstractDataGrid{D,T}

Stores the components of a wavefunction constructed from a planewave basis. Usually, the coefficient
data type `T` will be a `ComplexF32`, as in DFT calculations, double precision convergence of the
density will correspond to single-precision converegnce of the wavefunction.

Internally, coefficients are stored in an `Array{4,T}`. Indexing is then manually implemented, with
a `D`-dimensional `CartesianIndex` used for accessing each coefficient associated with a G-vector.
`PlanewaveWavefunction` instances are mutable, with `getindex()` and `setindex!()` defined for them,
but they are not resizable, and the backing array should not be resized.
"""
struct PlanewaveWavefunction{D,T} <: AbstractDataGrid{D,T}
    basis::ReciprocalBasis{D}
    spins::Vector{SVector{D,Float64}}
    kpoints::KPointList{D}
    energies::Array{Float64,3}
    occupancies::Array{Float64,3}
    grange::NTuple{D,UnitRange{Int}}
    data::Array{T,4}
    function PlanewaveWavefunction(
        basis::AbstractBasis{D},
        spins::AbstractVector{<:StaticVector{D,<:Real}},
        kpoints::AbstractKPointSet{D},
        energies::AbstractArray{<:Real,3},
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

"""
    PlanewaveWavefunction(
        [T::Type{<:Number} = Complex{Float32}],
        b::AbstractBasis{D},
        spins::AbstractVector{<:StaticVector{D,<:Real}},
        kpoints::AbstractKPointSet{D},
        bands::Integer
        grange::NTuple{D,<:AbstractUnitRange{<:Integer}},
    ) -> PlanewaveWavefunction{D,T}

Creates a new empty wavefunction with basis `b`, spin directions `spins`, a k-point list `kpoints`,
a number of bands `bands`, and a range of allowed G-vectors `grange`.
"""
function PlanewaveWavefunction(
    T::Type,
    basis::AbstractBasis{D},
    spins::AbstractVector{<:StaticVector{D,<:Real}},
    kpoints::AbstractKPointSet{D},
    bands::Integer,
    grange::NTuple{D,<:AbstractUnitRange{<:Integer}},
) where D
    energies = zeros(Float64, bands, length(kpoints), length(spins))
    occupancies = deepcopy(energies)
    data = zeros(T, prod(length.(grange)), bands, length(kpoints), length(spins))
    return PlanewaveWavefunction(basis, spins, kpoints, energies, occupancies, grange, data)
end

function PlanewaveWavefunction(
    basis::AbstractBasis{D},
    spins::AbstractVector{<:StaticVector{D,<:Real}},
    kpoints::AbstractKPointSet{D},
    bands::Integer,
    grange::NTuple{D,<:AbstractUnitRange{<:Integer}},
) where D
    return PlanewaveWavefunction(Complex{Float32}, basis, spins, kpoints, bands, grange)
end

Base.size(wf::PlanewaveWavefunction) = (reverse(size(wf.energies))..., length.(wf.grange)...)
Base.axes(wf::PlanewaveWavefunction) = (reverse(axes(wf.energies))..., wf.grange...)

PlanewaveIndices(wf::PlanewaveWavefunction) = PlanewaveIndices(axes(wf)...)

nspin(wf::PlanewaveWavefunction) = length(wf.spins)
nkpt(wf::PlanewaveWavefunction) = length(wf.kpoints)
nband(wf::PlanewaveWavefunction) = size(wf.energies, 1)

# Override some of the more generic AbstractDataGrid methods
Base.getindex(wf::PlanewaveWavefunction, i...) = throw(MethodError(getindex, (wf, i...)))
Base.setindex!(wf::PlanewaveWavefunction, i...) = throw(MethodError(setindex!, (wf, i...)))

function Base.getindex(wf::PlanewaveWavefunction{D}, i::PlanewaveIndex{D}) where D
    l = LinearIndices(wf.grange)[CartesianIndex(mod.(Tuple(i.g), length.(wf.grange)) .+ 1)]
    return wf.data[l, i.band, i.kpoint, i.spin]
end

# Broken due to the lack of a k-point data structure
function Base.getindex(wf::PlanewaveWavefunction, spin, kpt=:, band=:)
    return PlanewaveWavefunction(
        wf.basis,
        wf.spins[spin],
        wf.kpoints[kpt],
        wf.energies[band, kpt, spin],
        wf.occupancies[band, kpt, spin],
        wf.grange,
        wf.data[:, band, kpt, spin]
    )
end

function Base.getindex(
    wf::PlanewaveWavefunction{D},
    spin::Integer,
    kpt::Integer,
    band::Integer,
    g::Vararg{<:Integer,D}
) where D
    return wf[PlanewaveIndex(spin, kpt, band, g)]
end

function Base.setindex!(wf::PlanewaveWavefunction{D}, x, i::PlanewaveIndex{D}) where D
    l = LinearIndices(wf.grange)[CartesianIndex(mod.(Tuple(i.g), length.(wf.grange)) .+ 1)]
    wf.data[l, i.band, i.kpoint, i.spin] = x
end

function Base.setindex!(
    wf::PlanewaveWavefunction{D},
    x,
    spin::Integer,
    kpt::Integer,
    band::Integer,
    g::Vararg{<:Integer,D}
) where D
    wf[PlanewaveIndex(spin, kpt, band, g)] = x
end
