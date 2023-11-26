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
