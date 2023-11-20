"""
    Multiplicity{M} <: AbstractUnitRange{Rational{Int}}

Represents a valid range of spin states corresponding with multiplicity `M`. These index as if they 
are `UnitRange{Rational{Int}}` objects of the form `-(M - 1)//2:(M - 1)//2`, and can be converted to
`UnitRange` objects with the `UnitRange` constructor.

These types are singleton types with a numeric `Int` parameter, allowing for its use as a type
parameter in other types which require information on the spin multiplicity or allowed spin states.

# Examples
```
julia> UnitRange(SpinRange{3}())
1:1

julia> SpinRange{4}() == -3//2:3//2
true
```
"""
struct Multiplicity{M} <: AbstractUnitRange{Rational{Int}}
    Multiplicity{M}() where M = (@assert M > 0 "M must be a positive integer."; new{Int(M)}())
end

"""
    Multiplicity(M)

Constructs `Multiplicity{M}()`. This function is not type-stable unless the value of `M` is known at
compile time, similar to `Val(M)`.
"""
Multiplicity(M) = Multiplicity{Int(M)}()

Base.size(::Multiplicity{M}) where M = (M,)
Base.first(::Multiplicity{M}) where M = -(M-1)//2
Base.last(::Multiplicity{M}) where M = (M-1)//2

@inbounds function Base.getindex(s::Multiplicity, i::Integer)
    return i in eachindex(s) ? first(s) + (i - 1) : throw(BoundsError(s, i))
end

Base.UnitRange(s::Multiplicity) = first(s):last(s)
Base.show(io::IO, s::Multiplicity) = print(io, typeof(s), "()")

#---Spin bivectors---------------------------------------------------------------------------------#
"""
    Electrum._is_skew_symmetric(m::AbstractMatrix)

Returns `true` if a matrix is skew-symmetric (equivalently, antisymmetric).
"""
function _is_skew_symmetric(m::AbstractMatrix)
    size(m, 1) == size(m, 2) || return false
    for a in axes(m, 1)
        for b in (a+1):lastindex(m, 1)
            m[a,b] == -m[b,a] || return false
        end
    end
    return true
end

"""
    SpinBivector{D,T}

A skew-symmetric matrix representing a spin [bivector](https://en.wikipedia.org/wiki/Bivector) in 
`D` dimensions. No automatic normalization has been implemented for this type.

# Why the spin axis is a bivector

The identification of a spin direction with a vector is only possible in 3D. This is because the
axis of rotation is defined by the [Hodge dual](https://en.wikipedia.org/wiki/Hodge_star_operator)
of the rotation plane. The Hodge dual is a linear map relating `k`-dimensional objects in 
`D`-dimensional space to `D-k`-dimensional objects in the same space. 

A bivector is a 2-dimensional object, and in three dimensions, its Hodge dual is a vector. This is
not the case in any other number of dimensions. To allow for generality in describing spin, this
representation is used instead.
"""
struct SpinBivector{D,T} <: StaticMatrix{D,D,T}
    data::SVector{D,SVector{D,T}}
    function SpinBivector(m::StaticMatrix{D,D,T}) where {D,T}
        @assert _is_skew_symmetric(m) "Input matrix is not skew-symmetric."
        return new{D,T}(SVector{D}(eachcol(m)))
    end
end

SpinBivector{D}(m::AbstractMatrix) where D = SpinBivector(SMatrix{D,D}(m))
SpinBivector{D,T}(m::AbstractMatrix) where {D,T} = SpinBivector(SMatrix{D,D,T}(m))

function Base.getproperty(b::SpinBivector, s::Symbol)
    s === :matrix && return hcat(getfield(b, :data)...)
    return getfield(b, s)
end

Base.propertynames(::SpinBivector; private = false) = private ? (:data, :matrix) : (:matrix,)

Base.getindex(b::SpinBivector, i...) = getindex(b.matrix)
# Really only for resolving method ambiguities
Base.getindex(b::SpinBivector, i::Int) = getindex(b.matrix, i)
Base.getindex(b::SpinBivector, i::Int...) = getindex(b.matrix, i...)

Base.:(==)(a::SpinBivector, b::SpinBivector) = (a.matrix == b.matrix)
DataSpace(::Type{<:SpinBivector{D}}) where D = ByRealSpace{D}()
# Required for StaticArray subtypes
Tuple(b::SpinBivector) = Tuple(b.matrix)

# Constructors for taking wedge products implicitly
"""
    SpinBivector(u::StaticVector{D}, v::StaticVector{D})
    SpinBivector{D}(u::AbstractVector, v::AbstractVector)
    SpinBivector{D,T}(u::AbstractVector, v::AbstractVector)

Constructs a spin bivector from the wedge products of vectors `u` and `v`, which define a bivector.

The first constructor may be used if only one argument is a `StaticVector{D}` (the other will 
automatically be converted to the correct size).
"""
SpinBivector(u::StaticVector{D}, v::StaticVector{D}) where D = SpinBivector(u*v' - v*u')
# Convert one argument to a StaticVector if needed
SpinBivector(u::StaticVector{D}, v::AbstractVector) where D = SpinBivector(u, SVector{D}(v))
SpinBivector(u::AbstractVector, v::StaticVector{D}) where D = SpinBivector(SVector{D}(u), v)

function SpinBivector{D}(u::AbstractVector, v::AbstractVector) where D
    return SpinBivector(SVector{D}(u), SVector{D}(v))
end

function SpinBivector{D,T}(u::AbstractVector, v::AbstractVector) where {D,T}
    return SpinBivector(SVector{D,T}(u), SVector{D,T}(v))
end
