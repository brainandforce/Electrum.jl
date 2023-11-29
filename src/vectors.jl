#---Shift vectors----------------------------------------------------------------------------------#
"""
    ShiftVector{S<:BySpace,D,T} <: StaticVector{D,T}

A vector in fractional coordinates representing a shift of a lattice or lattice dataset from the
origin. This wraps a `SVector{D,T}` with an optional weight parameter of type `T` that may be useful
when working with symmetrical structures or k-points in the irreducible Brillouin zone. If it is not
explicitly set, it defaults to 1.

If constructors do not explicitly reference an element type, the element type is automatically
inferred by promoting the types of the arguments.

# Type aliases

`ShiftVector` is a general way of working with vectors which shift a lattice or data within it, and
for this reason we define an alias for representing k-points:

    const KPoint = ShiftVector{ByReciprocalSpace}
"""
struct ShiftVector{S<:BySpace,D,T<:Real} <: StaticVector{D,T}
    vector::SVector{D,T}
    weight::T
    function ShiftVector{S,D,T}(vector::StaticVector, weight::Real = oneunit(T)) where {S,D,T}
        return new(vector, weight)
    end
end

const KPoint = ShiftVector{ByReciprocalSpace}
@doc (@doc ShiftVector) KPoint

# Needed to resolve method ambiguities
ShiftVector{S,D,T}(::StaticArray, ::Real = 1) where {S,D,T} = error("Argument must be a vector.")
ShiftVector{S,D}(::StaticArray, ::Real = 1) where {S,D} = error("Argument must be a vector.")
ShiftVector{S}(::StaticArray, ::Real = 1) where S = error("Argument must be a vector.")

function ShiftVector{S,D}(vector::StaticVector, weight::Real = true) where {S,D}
    T = promote_type(eltype(vector), typeof(weight))
    return ShiftVector{S,D,T}(vector, weight)
end

function ShiftVector{S}(vector::StaticVector{D}, weight::Real = true) where {S,D}
    T = promote_type(eltype(vector), typeof(weight))
    return ShiftVector{S,D,T}(vector, weight)
end

function ShiftVector{S,D,T}(vector::AbstractVector, weight::Real = true) where {S,D,T}
    return ShiftVector{S,D,T}(SVector{D}(vector), weight)
end

function ShiftVector{S,D}(vector::AbstractVector, weight::Real = true) where {S,D}
    T = promote_type(eltype(vector), typeof(weight))
    return ShiftVector{S,D,T}(SVector{D}(vector), weight)
end

ShiftVector{S}(coord::Real...; weight::Real = 1) where S = ShiftVector{S}(SVector(coord), weight)

Base.hash(s::ShiftVector, h::UInt) = hash(s.vector, hash(s.weight, h))

function Base.:(==)(u::ShiftVector{S1}, v::ShiftVector{S2}) where {S1,S2}
    return (S1 === S2 && u.vector == v.vector && u.weight == v.weight)
end

Base.IndexStyle(::Type{<:ShiftVector}) = IndexLinear()
Base.getindex(s::ShiftVector, i::Int) = s.vector[i]

Base.Tuple(s::ShiftVector) = Tuple(s.vector)

Base.zero(::Type{ShiftVector{S,D}}) where {S,D} = ShiftVector{S}(zero(SVector{D,Bool}))
Base.zero(::Type{ShiftVector{S,D,T}}) where {S,D,T} = ShiftVector{S}(zero(SVector{D,T}))

"""
    weight(k::ShiftVector{S,D,T}) -> T

Returns the weight associated with a `ShiftVector`.
"""
weight(s::ShiftVector) = s.weight

# TODO: can we implement a remainder that excludes -0.5?
"""
    truncate(s::ShiftVector) -> ShiftVector

Moves a `ShiftVector` so that its values lie within the range [-1/2, 1/2]. The weight is preserved.
"""
Base.truncate(s::ShiftVector) = (typeof(s))(rem.(s.vector, 1, RoundNearest), s.weight)

DataSpace(::Type{<:ShiftVector{S,D}}) where {S,D} = S{D}()
ByCoordinate(::Type{<:ShiftVector{S,D}}) where {S,D} = ByFractionalCoordinate{D}()

Base.summary(io::IO, s::ShiftVector) = print(io, typeof(s), " with weight ", s.weight)
Base.show(io::IO, s::ShiftVector) = print(io, typeof(s), '(', s.vector, ", ", s.weight, ')')

function Base.show(io::IO, k::KPoint)
    print(io, KPoint, '(')
    join(io, k.vector, ", ")
    print(io, ", weight = ", weight(k), ')')
end
