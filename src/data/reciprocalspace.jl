"""
    KPointGrid{D} <: AbstractKPoints{D}

Contains a grid used to generate k-points during a calculation.

The grid itself is given as an `SMatrix{D,D,Int}`, and can be interpreted as a set of `D` vectors 
given in terms of the primitive basis. These vectors can alternatively be used to construct a 
supercell.

The shift of the k-point mesh off Γ is given as an `SVector{D,Float64}`.
"""
struct KPointGrid{D} <: AbstractKPoints{D}
    grid::SMatrix{D,D,Int}
    orig::SVector{D,Float64}
    function KPointGrid{D}(grid::AbstractMatrix{<:Integer}, orig::AbstractVector{<:Real}) where D
        # only allow positive values in the grid matrix
        @assert all(x -> x >= 0, grid) "negative values are disallowed in the grid matrix"
        # Keep shift inside the Brillouin zone
        orig = orig -  round.(orig)
        return new(grid, orig)
    end
end

"""
    KPointList{D} <: AbstractKPoints{D}

Contains a list of k-points and their associated weights. This is useful when describing an ordered
list of k-points that are not associated with a grid - for instance, the k-points used in band 
structure calculations. In the future, it will be possible to convert a `KPointGrid` to a
`KPointList`.

The weights of k-points are used to compensate for their placement on sites with point symmetry. If
no weights are provided, they are all assumed to be equal. Any explicit input of k-points weights
will be normalized such that their sum is 1, following the convention used by abinit.
"""
struct KPointList{D} <: AbstractKPoints{D}
    points::Vector{SVector{D,Float64}}
    weights::Vector{Float64}
    function KPointList(
        points::AbstractVector{<:SVector{D,<:Real}},
        weights::AbstractVector{<:Real}
    ) where D
        @assert length(points) == length(weights) "Number of k-points and weights do not match."
        return new{D}(points, weights / sum(weights))
    end
end

function KPointList{D}(
    points::AbstractVector{<:AbstractVector{<:Real}},
    weights::AbstractVector{<:Real}
) where D
    @assert length.(points) == D "k-points have the wrong dimensionality"
    svpoints = [SVector{D,Float64}(v) for v in points]
    return KPointList(svpoints, weights)
end

function KPointList{D}(points::AbstractVector{<:AbstractVector{<:Real}}) where D
    return KPointList{D}(points, ones(length(points)))
end

function KPointList(points::AbstractVector{<:SVector{D,<:Real}}) where D
    return KPointList(points, ones(length(points)))
end

# Get the k-point and its associated weight as a NamedTuple
Base.getindex(k::KPointList, i) = (kpt=k.points[i], weight=k.weights[i])

function Base.setindex!(k::KPointList, v::AbstractVector, i)
    k.points[i] = v
end

"""
    nkpt(k::KPointList{D}) -> Int

Gets the number of k-points in a `KPointList`.
"""
nkpt(k::KPointList) = length(k.points)
Base.length(k::KPointList) = length(k.points)

#= TODO: figure out how to get a k-point list from a grid
#  This would require also getting the correct k-point weights
function KPointList{D}(k::KPointGrid)

end
=#

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

Constructs a new `BandAtKPoint` from a vector containing tuples of energy and occupancy data
(in that order).
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
generate the data (as a `KPointList{D}`)and the band information at every k-point (as a 
`Vector{BandAtKPoint}`).
"""
struct BandStructure{D} <: AbstractReciprocalSpaceData{D}
    # k-points for which band data is defined
    kpts::KPointList{D}
    # Set of energy and occupancy data
    bands::Vector{BandAtKPoint}
    function BandStructure{D}(kpts::KPointList{D}, bands::Vector{BandAtKPoint}) where D
        @assert nkpt(kpts) == length(bands) "Incorrect number of k-points or band datasets."
        @assert _allsame(length(bands)) "Number of bands is inconsistent."
        return new(kpts, bands)
    end
end

"""
    BandStructure{D}(kpts::AbstractKPoints{D}, bands::AbstractVector{<:BandAtKPoint}) where D

Generates a new band structure from k-point information and a vector containing band information
at each k-point.
"""
function BandStructure{D}(kpts::AbstractKPoints{D}, bands::AbstractVector{<:BandAtKPoint}) where D
    return BandStructure{D}(kpts, bands)
end

# Get the pair of a k-point and associated band data
function Base.getindex(b::BandStructure{D}, inds...) where D
    return (b.kpts[inds...], b.bands[inds...])
end

nkpt(b::BandStructure{D}) where D = nkpt(b.kpts)
nband(b::BandStructure{D}) where D = nband(b.bands[1])

"""
    FatBands{D} <: AbstractReciprocalSpaceData{D}

Stores information relevant to plotting fatbands.

- FatBands.bands: matrix of energies at each [kpt, band].
- FatBands.projband: array of lm-decomposed band structure. [orbital, ion, band, kpt].
- FatBands.cband: array of complex-valued contributions to band structure.
"""
struct FatBands{D} <: AbstractReciprocalSpaceData{D}
    bands::Matrix{Float64}
    projband::Array{Float64,4}
    cband::Array{Complex{Float64},4}
end

"""
    HKLData{D,T} <: AbstractReciprocalSpaceData{D}

Stores information associated with specific sets of reciprocal lattice vectors. Data can be
accessed and modified using regular indexing, where indices may be negative.

Internally, the data is stored such that the zero frequency components are at the first indices
along that dimension. 
"""
struct HKLData{D,T} <: AbstractHKL{D,T}
    basis::ReciprocalBasis{D}
    data::Array{T,D}
    function HKLData(
        basis::AbstractBasis{D},
        data::AbstractArray{T,D}
    ) where {D,T}
        return new{D,T}(basis, data)
    end
end

function HKLData(
    data::AbstractArray{T,D},
    bounds::AbstractVector{<:AbstractRange{<:Integer}}
) where {D,T}
    return HKLData(zero(BasisVectors{D}), data, bounds)
end

"""
    basis(hkl::HKLData)

Returns the real-space basis associated with an `HKLData` object.
"""
basis(hkl::HKLData) = hkl.basis

"""
    grid(hkl::HKLData{D,T}) -> Array{T,D}

Returns the array that contains the reciprocal space data. Note that this causes information
about the index offset to be lost!
"""
grid(hkl::HKLData) = hkl.data

Base.size(g::HKLData) = size(g.data)
Base.length(g::HKLData) = length(g.data)
Base.iterate(g::HKLData, i::Integer = 1) = iterate(g.grid, i)

# HKLData now supports indexing by Miller index
Base.getindex(g::HKLData, inds...) = getindex(g.grid, (mod.(inds, size(g)) .+ 1)...)

function Base.setindex!(g::HKLData, x, inds...)
    i = mod.(inds, size(g)) .+ 1
    g.data[i...] = x
end

Base.abs(hkl::HKLData) = HKLData(basis(hkl), abs.(grid(hkl)))
Base.abs2(hkl::HKLData) = HKLData(basis(hkl), abs2.(grid(hkl)))

"""
    voxelsize(g::HKLData)

Gets the size of a voxel asssociated with the `RealSpaceDataGrid` that would be generated by 
performing an inverse Fourier transform on the `HKLData`.
"""
voxelsize(g::HKLData) = volume(RealBasis(basis(g))) / length(g)

"""
    HKLDict{D,T}

An alternative to `HKLData` uses a dictionary instead of an array as a backing field.

This is a more space-efficient alternative to `HKLData` in the case of reciprocal space data with
a large number of zero components. For wavefunction data, which is often specified to some energy
cutoff that corresponds to a distance in reciprocal space, there are many zero valued elements to
the array. Unspecified elements in an `HKLDict` are assumed to be zero.
"""
struct HKLDict{D,T} <: AbstractHKL{D,T}
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

"""
    vectors(hkl::HKLDict)

Returns the set of vectors in an `HKLDict` for which values have been defined.
"""
function vectors(hkl::HKLDict)
    return keys(hkl.dict)
end

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
    ReciprocalWavefunction{D,T<:Real} <: AbstractReciprocalSpaceData{D}

Contains a wavefunction stored by k-points and bands in a planewave basis. Used to store data in
VASP WAVECAR files. Each k-point is expected to have the same number of bands.

Every band has associated data containing coefficients of the constituent planewaves stored in a 
`HKLData{D,Complex{T}}`. Unlike most data structures provided by this package, the type of
complex number used does not default to `Float64`: wavefunction data is often supplied as a 
`Complex{Float32}` since wavefunctions usually only converge to single precision, and `Float64`
storage would waste space.
"""
struct ReciprocalWavefunction{D,T<:Real} <: AbstractReciprocalSpaceData{D}
    # Reciprocal lattice on which the k-points are defined
    rlatt::ReciprocalBasis{D}
    # k-points used to construct the wavefunction
    kpts::KPointList{D}
    # Planewave coefficients: an Array{HKLData,3} (size nspin*nkpt*maxnband)
    waves::Array{HKLData{D,Complex{T}},3}
    # Energies and occupancies, Array{Float64,3} with the same size as above
    energies::Array{Float64,3}
    occupancies::Array{Float64,3}
    function ReciprocalWavefunction(
        rlatt::AbstractBasis{D},
        kpts::AbstractKPoints{D},
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
    kpts::AbstractKPoints{D},
    waves::AbstractArray{HKLData{D,Complex{T}},3}
) where {D,T<:Real}
    # Construct zero matrix
    z = zeros(Float64, size(waves))
    return ReciprocalWavefunction(rlatt, kpts, waves, z, z)
end

function ReciprocalWavefunction(
    latt::AbstractLattice{D},
    kpts::AbstractKPoints{D},
    waves::AbstractArray{HKLData{D,Complex{T}},3},
    energies::AbstractArray{<:Real,3},
    occupancies::AbstractArray{<:Real,3};
    primitive::Bool = true
) where {D,T<:Real}
    M = primitive ? prim(latt) : conv(latt)
    return ReciprocalWavefunction(M, kpts, waves, energies, occupancies)
end

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
