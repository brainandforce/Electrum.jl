"""
    KPoint{D} <: DenseVector{D,Float64}

Stores a k-point with an associated weight that corresponds to the number of symmetry-equivalent
k-points, stored as an integer.
"""
struct KPoint{D} <: DenseVector{Float64}
    point::SVector{D,Float64}
    weight::Int
    KPoint(pt::StaticVector{D,<:Real}, wt::Integer = 1) where D = new{D}(pt .- round.(pt), wt)
end

KPoint(x::Real...; weight::Integer = 1) = KPoint(SVector(x), weight)
KPoint{D}(pt::AbstractVector{<:Real}, wt::Integer = 1) where D = KPoint(SVector{D,Float64}(pt), wt)

Base.hash(k::KPoint, h::UInt) = hash(hash(k.point, hash(k.weight)), h)
Base.:(==)(k1::KPoint, k2::KPoint) = k1.point == k2.point && k1.weight == k2.weight

Base.length(::Type{KPoint{D}}) where D = D
Base.size(::Type{KPoint{D}}) where D = (D,)
Base.size(k::KPoint) = size(k.point)
Base.axes(::Type{KPoint{D}}) where D = (Base.OneTo(D),)
Base.axes(k::KPoint) = axes(k.point)

Base.IndexStyle(::Type{<:KPoint}) = IndexLinear()
Base.getindex(k::KPoint, i) = k.point[i]

Base.convert(T::Type{<:StaticVector}, k::KPoint) = convert(T, k.point)::T
Base.convert(T::Type{<:Vector}, k::KPoint) = convert(T, k.point)::T
Base.convert(T::Type{<:KPoint}, v::AbstractVector{<:Real}) = T(v, weight = 1)
Base.convert(T::Type{<:KPoint}, k::KPoint) = k::T

Base.zero(::Type{KPoint{D}}) where D = KPoint(zero(SVector{D,Float64}))

"""
    weight(k::KPoint) -> Int

Returns the weight associated with a k-point.
"""
weight(k::KPoint) = k.weight

#---Generated lists of k-points--------------------------------------------------------------------#

"""
    KPointMesh{D} <: AbstractVector{KPoint{D}}

Contains a list of k-points associated with a matrix describing the mesh that was used to generate
the points, and its shift off the Γ point (origin). If the mesh used to generate the points is
unknown, it will be set to the zero matrix of dimension `D`.

A `KPointMesh` can be indexed as if it were an ordinary `Vector{KPoint{D}}`.
"""
struct KPointMesh{D} <: AbstractVector{KPoint{D}}
    points::Vector{KPoint{D}}
    grid::SMatrix{D,D,Int}
    shift::SVector{D,Float64}
    function KPointMesh(
        points::AbstractVector{KPoint{D}}, 
        grid::StaticMatrix{D,D,<:Integer} = zeros(SMatrix{D,D,Int}),
        shift::StaticVector{D,<:Real} = zeros(SVector{D,Float64})
    ) where D
        return new{D}(points, grid, shift)
    end
end

function KPointMesh(
    points::AbstractVector,
    grid::StaticMatrix{D,D,<:Integer} = zeros(SMatrix{D,D,Int}),
    shift::StaticVector{D,<:Real} = zeros(SVector{D,Float64})
) where D
    return KPointMesh(KPoint{D}.(points), grid, shift)
end

function KPointMesh{D}(
    points::AbstractVector,
    grid::AbstractMatrix{<:Integer} = zeros(SMatrix{D,D,Int}),
    shift::AbstractVector{<:Real} = zeros(SVector{D,Float64})
) where D
    return KPointMesh(points, SMatrix{D,D,Int}(grid), SVector{D,Float64}(shift))
end

function Base.:(==)(k1::KPointMesh, k2::KPointMesh)
    return k1.points == k2.points && k1.grid == k2.grid && k1.shift == k2.shift
end

Base.size(k::KPointMesh) = size(k.points)
Base.axes(k::KPointMesh) = axes(k.points)
Base.IndexStyle(::Type{<:KPointMesh}) = IndexLinear()
Base.getindex(k::KPointMesh, i) = k.points[i]
Base.iterate(k::KPointMesh, i::Integer = 1) = iterate(k.points, i)
Base.convert(T::Type{Vector{<:KPoint}}, k::KPointMesh) = k.points::T

"""
    BandAtKPoint

Stores information about a band's energy and its occupancy at a specific k-point.
"""
struct BandAtKPoint
    # Energies
    e::Vector{Float64}
    # Occupancy
    occ::Vector{Float64}
    function BandAtKPoint(e::AbstractVector{<:Real}, occ::AbstractVector{<:Real})
        @assert length(e) == length(occ) "Size of energy and occupancy arrays do not match."
        return new(e, occ)
    end
end

"""
    BandAtKPoint(eocc::AbstractVector{NTuple{2,<:Real}})

Constructs a new `BandAtKPoint` from a vector containing tuples of energy and occupancy data (in 
that order).
"""
function BandAtKPoint(eocc::AbstractVector{<:NTuple{2,<:Real}})
    return BandAtKPoint([x[1] for x in eocc], [x[2] for x in eocc])
end

# Access any pair of energy and occupancy with indexing
function Base.getindex(b::BandAtKPoint, inds...)
    return (b.e[inds...], b.occ[inds...])
end

"""
    nband(b::BandAtKPoint) -> Int

Returns the number of bands associated with a k-point.
"""
nband(b::BandAtKPoint) = length(b.e)

"""
    BandStructure{D}

Stores information about an electronic band structure, including the list of k-points used to
generate the data (as a `KPointMesh{D}`)and the band information at every k-point (as a
`Vector{BandAtKPoint}`).
"""
struct BandStructure{D}
    # k-points for which band data is defined
    kpts::Vector{KPoint{D}}
    # Set of energy and occupancy data
    bands::Vector{BandAtKPoint}
    function BandStructure{D}(kpts::AbstractVector{KPoint{D}}, bands::Vector{BandAtKPoint}) where D
        @assert nkpt(kpts) == length(bands) "Incorrect number of k-points or band datasets."
        @assert _allsame(length(bands)) "Number of bands is inconsistent."
        return new(kpts, bands)
    end
end

"""
    BandStructure{D}(kpts::AbstractKPoints{D}, bands::AbstractVector{<:BandAtKPoint}) where D

Generates a new band structure from k-point information and a vector containing band information at
each k-point.
"""
function BandStructure{D}(kpts::AbstractKPointSet{D}, bands::AbstractVector{<:BandAtKPoint}) where D
    return BandStructure{D}(kpts, bands)
end

# Get the pair of a k-point and associated band data
function Base.getindex(b::BandStructure{D}, inds...) where D
    return (b.kpts[inds...], b.bands[inds...])
end

nkpt(b::BandStructure{D}) where D = nkpt(b.kpts)
nband(b::BandStructure{D}) where D = nband(b.bands[1])

"""
    FatBands{D}

Stores information relevant to plotting fatbands.

- FatBands.bands: matrix of energies at each [kpt, band].
- FatBands.projband: array of lm-decomposed band structure. [orbital, ion, band, kpt].
- FatBands.cband: array of complex-valued contributions to band structure.
"""
struct FatBands{D}
    bands::Matrix{Float64}
    projband::Array{Float64,4}
    cband::Array{Complex{Float64},4}
end

"""
    HKLData{D,T} <: AbstractDataGrid{D,T}

Stores information associated with a reciprocal space basis. Data can be accessed and modified by
using the G-vectors as indices. Associated k-point data is also provided; if no k-point is supplied
it is assumed to be the gamma point (`zero(SVector{D,Float64})`)

Internally, the data is stored such that the zero frequency components are at the first indices
along that dimension. The data at G-vector `[0, 0, 0]` is stored in the backing array's `[1, 1, 1]`
index, and the rest of the indices correspond to reciprocal space points using the FFT convention.
"""
struct HKLData{D,T} <: AbstractDataGrid{D,T}
    basis::ReciprocalBasis{D}
    data::Array{T,D}
    kpt::KPoint{D}
    function HKLData(
        basis::AbstractBasis{D},
        data::AbstractArray{T,D},
        kpt::AbstractVector{<:Real} = zero(KPoint{D})
    ) where {D,T}
        # TODO: Do we want to perform circular shifts of the data?
        # How does this work with array copying?
        # Move k-point so the indices are in the range (-0.5, 0.5]
        return new{D,T}(basis, data, kpt)
    end
end

function grid_specific_check(g::HKLData...)
    any(h -> !isapprox(first(g).kpt, h.kpt), g) && error("k-points do not match.")
    return nothing
end

function Base.zeros(
    ::Type{HKLData{D,T}},
    basis::AbstractBasis{D},
    ranges::Vararg{AbstractUnitRange{<:Integer},D},
) where {D,T}
    return HKLData(ReciprocalBasis(basis), zeros(T, length.(ranges)))
end

Base.abs(hkl::HKLData) = HKLData(basis(hkl), abs.(g.data))
Base.abs2(hkl::HKLData) = HKLData(basis(hkl), abs2.(g.data))

"""
    voxelsize(g::HKLData)

Gets the size of a voxel asssociated with the `RealSpaceDataGrid` that would be generated by 
performing an inverse Fourier transform on the `HKLData`.
"""
voxelsize(g::HKLData) = volume(RealBasis(g)) / length(g)

function Base.isapprox(g1::HKLData, g2::HKLData; kwargs...)
    @assert basis(g1) === basis(g2) "Grid basis vectors for each grid are not identical."
    @assert size(g.data) === size(g.data) "Grid sizes are different."
    return isapprox(g.data, g.data, kwargs...)
end

"""
    HKLDict{D,T}

An alternative to `HKLData` uses a dictionary instead of an array as a backing field.

This is a more space-efficient alternative to `HKLData` in the case of reciprocal space data with a
large number of zero components. For wavefunction data, which is often specified to some energy
cutoff that corresponds to a distance in reciprocal space, there are many zero valued elements to
the array. Unspecified elements in an `HKLDict` are assumed to be zero.
"""
struct HKLDict{D,T} <: AbstractDataGrid{D,T}
    dict::Dict{SVector{D,Int},T}
end

Base.has_offset_axes(hkl::HKLDict) = true

function Base.getindex(hkl::HKLDict{D,T}, inds...) where {D,T}
    v = SVector{D,Int}(inds...)
    if haskey(hkl.dict, v)
        return hkl.dict[v]
    else
        # Return a zero element of some kind by default
        return zero(T)
    end
end

function Base.setindex!(hkl::HKLDict{D,T}, value::T, inds...) where {D,T}
    hkl.dict[SVector{D,Int}(inds...)] = value
end

Base.keys(hkl::HKLDict) = keys(hkl.dict)

Base.iterate(hkl::HKLDict) = iterate(hkl.dict)
Base.iterate(hkl::HKLDict, i) = iterate(hkl.dict, i)

function HKLDict(hkl::HKLData{D,T}) where {D,T<:Union{<:Number,<:AbstractArray{Number}}}
    dict = Dict{SVector{D,Int},T}()
    # Get the offset for the indices
    offset = minimum.(hkl.bounds) .- 1
    # Iterate through the matrix and get its coordinates
    for ind in CartesianIndices(hkl.data)
        # Only add nonzero elements
        if hkl.data[ind] != zero(T)
            dict[Tuple(ind) .+ offset] = hkl.data[ind]
        end
    end
    return HKLDict(dict)
end

function HKLData(hkl::HKLDict{D,T}) where {D,T<:Union{<:Number,<:AbstractArray{Number}}}
    # Find the bounds
    bounds = MVector(UnitRange(extrema(v[n] for v in keys(hkl.dict)...)) for n in 1:D)
    data = zeros(T, length.(bounds)...)
    # Loop through the dictionary
    for (k,v) in hkl.dict
        data[k...] = v
    end
    return HKLData(data, bounds)
end

"""
    ReciprocalWavefunction{D,T<:Real}

Contains a wavefunction stored by k-points and bands in a planewave basis. Used to store data in
VASP WAVECAR files. Each k-point is expected to have the same number of bands.

Every band has associated data containing coefficients of the constituent planewaves stored in a
`HKLData{D,Complex{T}}`. Unlike most data structures provided by this package, the type of complex
number used does not default to `Float64`: wavefunction data is often supplied as a 
`Complex{Float32}` since wavefunctions usually only converge to single precision, and `Float64`
storage would waste space.

The energies and occupancies are also stored in fields with the corresponding names, and can be
accessed by spins, k-points, and bands, with indices in that order.
"""
struct ReciprocalWavefunction{D,T<:Real}
    # Reciprocal lattice on which the k-points are defined
    rlatt::ReciprocalBasis{D}
    # k-points used to construct the wavefunction
    kpts::KPointMesh{D}
    # Planewave coefficients: an Array{HKLData,3} (size nspin*nkpt*maxnband)
    waves::Array{HKLData{D,Complex{T}},3}
    # Energies and occupancies, Array{Float64,3} with the same size as above
    energies::Array{Float64,3}
    occupancies::Array{Float64,3}
    function ReciprocalWavefunction(
        rlatt::AbstractBasis{D},
        kpts::AbstractVector{KPoint{D}},
        waves::AbstractArray{HKLData{D,Complex{T}},3},
        energies::AbstractArray{<:Real,3},
        occupancies::AbstractArray{<:Real,3},
    ) where {D,T<:Real}
        @assert length(kpts) == size(waves, 2) string(
            "k-point list length inconsistent with number of wavefunction entries"
        )
        return new{D,T}(rlatt, kpts, waves, energies, occupancies)
    end
end

# When eneregies and occupancies are not specified
function ReciprocalWavefunction(   
    rlatt::AbstractBasis{D},
    kpts::AbstractVector{KPoint{D}},
    waves::AbstractArray{HKLData{D,Complex{T}},3}
) where {D,T<:Real}
    # Construct zero matrix
    z = zeros(Float64, size(waves))
    return ReciprocalWavefunction(rlatt, kpts, waves, z, z)
end

data_space(::Type{<:ReciprocalWavefunction{D}}) where D = ByReciprocalSpace{D}()

"""
    bounds(wf::ReciprocalWavefunction)

Gets the range of valid G-vectors in a `ReciprocalWavefunction`.
"""
function bounds(wf::ReciprocalWavefunction{D,T}) where {D,T}
    inds = CartesianIndices((0:0, 0:0, 0:0))
    # Loop through each HKLData
    for hkl in wf.waves
        i = CartesianIndices(hkl)
        # Skip this if the new indices are equal
        if i != inds
            # Create a tuple with every longer range
            inds = CartesianIndices(
                NTuple{D,UnitRange{Int}}(
                    length(a) >= length(b) ? a : b
                    for (a,b) in zip(i.indices, inds.indices)
                )
            )
        end
    end
    return inds
end

Base.size(wf::ReciprocalWavefunction) = size(wf.waves)
Base.length(wf::ReciprocalWavefunction) = length(wf.waves)

function Base.getindex(wf::ReciprocalWavefunction, inds...)
    return (
        coeffs = wf.waves[inds...],
        energies = wf.energies[inds...],
        occupancies = wf.occupancies[inds...]
    )
end

"""
    nspin(wf::ReciprocalWavefunction) -> Int

Returns the number of spins associated with a `ReciprocalWavefunction`.
"""
nspin(wf::ReciprocalWavefunction) = size(wf.waves, 1)

"""
    nkpt(wf::ReciprocalWavefunction) -> Int

Returns the number of k-points associated with a `ReciprocalWavefunction`.
"""
nkpt(wf::ReciprocalWavefunction) = size(wf.waves, 2)

"""
    nband(wf::ReciprocalWavefunction) -> Int

Returns the number of bands associated with a `ReciprocalWavefunction`. It is assumed that the
number of bands is the same for each k-point and spin.
"""
nband(wf::ReciprocalWavefunction) = size(wf.waves, 3)

basis(wf::ReciprocalWavefunction) = wf.rlatt

"""
    fermi(wf::ReciprocalWavefunction) -> Float64

Estimates the Fermi energy associated with a reciprocal space wavefunction using the energy and
occupancy data in the `ReciprocalWavefunction`.
"""
function fermi(wf::ReciprocalWavefunction)
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
