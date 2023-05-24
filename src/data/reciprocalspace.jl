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

Base.hash(k::KPoint, h::UInt) = hash(k.point, hash(k.weight, h))
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
Base.convert(T::Type{<:KPoint}, v::AbstractVector{<:Real}) = T(v, 1)
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
Base.setindex!(k::KPointMesh, x, i) = setindex!(k.points, x, i)

Base.convert(T::Type{Vector{<:KPoint}}, k::KPointMesh) = k.points::T
Base.convert(T::Type{<:KPointMesh}, v::AbstractVector{<:KPoint}) = KPointMesh(v)::T
Base.convert(T::Type{<:KPointMesh}, k::KPointMesh) = k::T

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
generate the data (as am `AbstractVector{KPoint{D}}`)and the band information at every k-point (as a
`Vector{BandAtKPoint}`).
"""
struct BandStructure{D}
    # k-points for which band data is defined
    kpts::Vector{KPoint{D}}
    # Set of energy and occupancy data
    bands::Vector{BandAtKPoint}
    function BandStructure(kpts::AbstractVector{KPoint{D}}, bands::Vector{BandAtKPoint}) where D
        @assert nkpt(kpts) == length(bands) "Incorrect number of k-points or band datasets."
        @assert _allsame(length(bands)) "Number of bands is inconsistent."
        return new{D}(kpts, bands)
    end
end

"""
    BandStructure(kpts::AbstractVector{KPoint{D}}, bands::AbstractVector{<:BandAtKPoint})

Generates a new band structure from k-point information and a vector containing band information at
each k-point.
"""
function BandStructure(
    kpts::AbstractVector{KPoint{D}},
    bands::AbstractVector{<:BandAtKPoint}
) where D
    return BandStructure(kpts, bands)
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
