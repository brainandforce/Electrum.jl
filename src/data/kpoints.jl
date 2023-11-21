"""
    KPoint{D,T<:Real} <: StaticVector{D,T}

Stores a k-point as reduced reciprocal space coordiantes with an associated weight that corresponds
to the number of symmetry-equivalent k-points, stored as an integer.
"""
struct KPoint{D,T<:Real} <: StaticVector{D,T}
    point::SVector{D,T}
    weight::Int
    KPoint{D,T}(pt::StaticVector, wt::Integer = 1) where {D,T} = new(pt, wt)
end

# Needed to resolve method ambiguities
KPoint{D,T}(::StaticArray, ::Integer = 1) where {D,T} = error("Argument must be a vector.")
KPoint{D}(::StaticArray, ::Integer = 1) where D = error("Argument must be a vector.")
KPoint(::StaticArray, ::Integer = 1) = error("Argument must be a vector.")

KPoint{D}(pt::StaticVector, wt::Integer = 1) where D = KPoint{D,eltype(pt)}(pt, wt)
KPoint(pt::StaticVector, wt::Integer = 1) = KPoint{length(pt),eltype(pt)}(pt, wt)

KPoint{D,T}(pt::AbstractVector, wt::Integer = 1) where {D,T} = KPoint(SVector{D,T}(pt), wt)
KPoint{D}(pt::AbstractVector, wt::Integer = 1) where D = KPoint(SVector{D}(pt), wt)

KPoint(pt::Real...; weight::Integer = 1) = KPoint(SVector(pt), weight)

Base.hash(k::KPoint, h::UInt) = hash(k.point, hash(k.weight, h))
Base.:(==)(k1::KPoint, k2::KPoint) = k1.point == k2.point && k1.weight == k2.weight

Base.IndexStyle(::Type{<:KPoint}) = IndexLinear()
Base.getindex(k::KPoint, i::Int) = k.point[i]

Tuple(k::KPoint) = Tuple(k.point)
# Base.convert(T::Type{<:AbstractVector}, k::KPoint) = convert(T, k.point)

Base.zero(::Type{KPoint{D}}) where D = KPoint(zero(SVector{D,Bool}))
Base.zero(::Type{KPoint{D,T}}) where {D,T} = KPoint(zero(SVector{D,T}))

"""
    weight(k::KPoint) -> Int

Returns the weight associated with a k-point.
"""
weight(k::KPoint) = k.weight

# TODO: can we implement a remainder that excludes -0.5?
"""
    truncate(k::KPoint) -> KPoint

Moves a k-point so that its values lie within the range [-1/2, 1/2]. The weight is preserved.
"""
Base.truncate(k::KPoint) = KPoint(rem.(k.point, 1, RoundNearest), k.weight)

#---Generated lists of k-points--------------------------------------------------------------------#
"""
    KPointMesh{D,T} <: AbstractVector{KPoint{D,T}}

Contains a list of k-points associated with a matrix describing the mesh that was used to generate
the points, and its shift off the Γ point (origin). If the mesh used to generate the points is
unknown, it will be set to the zero matrix of dimension `D`.

A `KPointMesh` can be indexed as if it were an ordinary `Vector{KPoint{D,T}}`.
"""
struct KPointMesh{D,T} <: AbstractVector{KPoint{D,T}}
    points::Vector{KPoint{D,T}}
    grid::SMatrix{D,D,Int}
    shift::SVector{D,T}
    function KPointMesh(
        points::AbstractArray{KPoint{D,T}}, 
        grid::StaticMatrix{D,D} = zeros(SMatrix{D,D,Int}),
        shift::StaticVector{D} = zeros(SVector{D,T})
    ) where {D,T}
        return new{D,promote_type(T, eltype(shift))}(points, grid, shift)
    end
end

function KPointMesh(
    points::AbstractArray,
    grid::StaticMatrix{D,D} = zeros(SMatrix{D,D,Int}),
    shift::StaticVector{D} = zeros(SVector{D,Bool})
) where D
    return KPointMesh(KPoint{D}.(points), grid, shift)
end

function KPointMesh{D}(
    points::AbstractArray,
    grid::AbstractMatrix = zeros(SMatrix{D,D,Int}),
    shift::AbstractVector = zeros(SVector{D,Bool})
) where D
    return KPointMesh(KPoint{D}.(points), SMatrix{D,D,Int}(grid), SVector{D,T}(shift))
end

function KPointMesh{D,T}(
    points::AbstractArray,
    grid::AbstractMatrix = zeros(SMatrix{D,D,Int}),
    shift::AbstractVector = zeros(SVector{D,T})
) where {D,T}
    return KPointMesh(KPoint{D,T}.(points), SMatrix{D,D,Int}(grid), SVector{D,T}(shift))
end

function Base.:(==)(k1::KPointMesh, k2::KPointMesh)
    return k1.points == k2.points && k1.grid == k2.grid && k1.shift == k2.shift
end

Base.size(k::KPointMesh) = size(k.points)
Base.axes(k::KPointMesh) = axes(k.points)

Base.IndexStyle(::Type{<:KPointMesh}) = IndexLinear()
Base.getindex(k::KPointMesh, i) = k.points[i]
Base.setindex!(k::KPointMesh, x, i) = setindex!(k.points, x, i)
Base.setindex!(k::KPointMesh, v::AbstractVector, i) = setindex!(k.points, eltype(k)(v), i)

Base.convert(T::Type{Vector{<:KPoint}}, k::KPointMesh) = k.points::T
Base.convert(T::Type{<:KPointMesh}, v::AbstractVector{<:KPoint}) = KPointMesh(v)::T
Base.convert(T::Type{<:KPointMesh}, k::KPointMesh) = k::T

"""
    nkpt(k::KPointMesh) -> Int

Counts the number of k-points (equivalent to `length(k)`). This function can be defined for custom
types that contain a `KPointMesh`.
"""
nkpt(k::KPointMesh) = length(k)

"""
    nkpt(x) -> Int

Counts the number of explicitly enumerated k-points associated with an object containing a
`KPointMesh`.

By default, this function returns `length(KPointMesh(x))`, so defining `KPointMesh(x::T)` for a
type `T` containing a `KPointMesh` will automatically define `nkpt(x::T)`.
"""
nkpt(x) = length(KPointMesh(x))
