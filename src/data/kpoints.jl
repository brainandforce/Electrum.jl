"""
    KPoint{D,T} <: DenseVector{D,T}

Stores a k-point with an associated weight that corresponds to the number of symmetry-equivalent
k-points, stored as an integer.
"""
struct KPoint{D,T} <: DenseVector{T}
    point::SVector{D,T}
    weight::Int
    KPoint{D,T}(pt::AbstractVector{<:Real}, wt::Integer = 1) where {D,T} = new(pt .- round(pt), wt)
end

KPoint{D}(pt::AbstractVector{T}, wt::Integer = 1) where {D,T} = KPoint{D,T}(pt, wt)
KPoint(pt::SVector{D,T}, wt::Integer = 1) where {D,T} = KPoint{D,T}(pt, wt)

KPoint(x::Real...; weight::Integer = 1) = KPoint(SVector(x), weight)

Base.hash(k::KPoint, h::UInt) = hash(k.point, hash(k.weight, h))
Base.:(==)(k1::KPoint, k2::KPoint) = k1.point == k2.point && k1.weight == k2.weight

Base.size(::Type{<:KPoint{D}}) where D = (D,)
Base.size(k::KPoint) = size(k.point)
Base.axes(::Type{<:KPoint{D}}) where D = (Base.OneTo(D),)
Base.axes(k::KPoint) = axes(k.point)

Base.IndexStyle(::Type{<:KPoint}) = IndexLinear()
Base.getindex(k::KPoint, i) = k.point[i]

Base.convert(T::Type{<:StaticVector}, k::KPoint) = convert(T, k.point)::T
Base.convert(T::Type{<:Vector}, k::KPoint) = convert(T, k.point)::T

Base.convert(T::Type{<:KPoint}, v::AbstractVector{<:Real}) = T(v, 1)
Base.convert(T::Type{<:KPoint}, k::KPoint) = k::T

Base.zero(::Type{KPoint{D,T}}) where {D,T} = KPoint(zero(SVector{D,T}))

"""
    weight(k::KPoint) -> Int

Returns the weight associated with a k-point.
"""
weight(k::KPoint) = k.weight

#---Generated lists of k-points--------------------------------------------------------------------#

"""
    KPointMesh{D,T} <: AbstractVector{KPoint{D,T}}

Contains a list of k-points associated with a matrix describing the mesh that was used to generate
the points, and its shift off the Γ point (origin). If the mesh used to generate the points is
unknown, it will be set to the zero matrix of dimension `D`.

A `KPointMesh` can be indexed as if it were an ordinary `Vector{KPoint{D}}`.
"""
struct KPointMesh{D,T} <: AbstractVector{KPoint{D,T}}
    points::Vector{KPoint{D,T}}
    grid::SMatrix{D,D,Int}
    shift::SVector{D,T}
    function KPointMesh{D,T}(
        points::AbstractVector{<:AbstractVector},
        grid::AbstractMatrix{<:Integer} = zero(SMatrix{3,3,Int}),
        shift::AbstractVector = zero(SVector{D,T})
    ) where {D,T}
        return new{D,T}(points, grid, shift)
    end
end

function KPointMesh{D}(
    points::AbstractVector{<:AbstractVector{T}},
    grid::AbstractMatrix{<:Integer} = zeros(SMatrix{D,D,Int}),
    shift::AbstractVector{<:Real} = zeros(SVector{D,Bool})
) where {D,T}
    return KPointMesh{D,promote_type(T, eltype(shift))}(points, grid, shift)
end

function KPointMesh(
    points::AbstractVector{<:AbstractVector{T}},
    grid::StaticMatrix{D,D,<:Integer} = zeros(SMatrix{D,D,Int}),
    shift::StaticVector{D,<:Real} = zeros(SVector{D,Bool})
) where {D,T}
    return KPointMesh{D,promote_type(T, eltype(shift))}(points, grid, shift)
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
    nkpt(k::KPointMesh) -> Int

Counts the number of k-points (equivalent to `length(k)`). This function can be defined for custom
types that contain a `KPointMesh`.
"""
nkpt(k::KPointMesh) = length(k)
