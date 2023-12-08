#---Coordinate vectors-----------------------------------------------------------------------------#
"""
    CoordinateVector{S<:BySpace,C<:ByCoordinate,D,T<:Real} <: StaticVector{D,T}

Represents a spatial coordinate with space given by trait `S` and coordinate trait `C`. Coordinate
values must be subtypes of `Real`.
"""
struct CoordinateVector{S<:BySpace,C<:ByCoordinate,D,T<:Real} <: StaticVector{D,T}
    vector::SVector{D,T}
    CoordinateVector{S,C}(v::StaticVector{D,T}) where {S,C,D,T} = new{S,C,D,T}(v)
end

BySpace(::Type{<:CoordinateVector{S}}) where S = S()
ByCoordinate(::Type{<:CoordinateVector{<:BySpace,C}}) where C = C()

"""
    RealCartesianCoordinate{D,T}
    (alias for CoordinateVector{Electrum.ByRealSpace,Electrum.ByCartesianCoordinate,D,T})

Represents a real space Cartesian coordinate. For more information, see [`CoordinateVector`](@ref).
"""
const RealCartesianCoordinate = CoordinateVector{ByRealSpace,ByCartesianCoordinate}

"""
    RealFractionalCoordinate{D,T}
    (alias for CoordinateVector{Electrum.ByRealSpace,Electrum.ByFractionalCoordinate,D,T})

Represents a real space fractional coordinate. For more information, see [`CoordinateVector`](@ref).
"""
const RealFractionalCoordinate = CoordinateVector{ByRealSpace,ByFractionalCoordinate}

"""
    ReciprocalCartesianCoordinate{D,T}
    (alias for CoordinateVector{Electrum.ByReciprocalSpace,Electrum.ByCartesianCoordinate,D,T})

Represents a reciprocal space Cartesian coordinate. For more information, see
[`CoordinateVector`](@ref).
"""
const ReciprocalCartesianCoordinate = CoordinateVector{ByReciprocalSpace,ByCartesianCoordinate}

"""
    ReciprocalFractionalCoordinate{D,T}
    (alias for CoordinateVector{Electrum.ByReciprocalSpace,Electrum.ByFractionalCoordinate,D,T})

Represents a reciprocal space fractional coordinate. For more information,
see [`CoordinateVector`](@ref).
"""
const ReciprocalFractionalCoordinate = CoordinateVector{ByReciprocalSpace,ByFractionalCoordinate}

array_not_flattened() = error(
    "Multidimensional array arguments are not automatically flattened to vectors with this " *
    "constructor. To do this, explicitly include the dimension type parameter."
)

CoordinateVector{S,C}(::StaticArray) where {S,C} = array_not_flattened()
CoordinateVector{S,C}(t::Tuple) where {S,C} = CoordinateVector{S,C}(SVector(t))
CoordinateVector{S,C}(x::Real...) where {S,C} = CoordinateVector{S,C}(SVector(x))

CoordinateVector{S,C,D}(t::Tuple) where {S,C,D} = CoordinateVector{S,C}(SVector(t))
CoordinateVector{S,C,D}(v::StaticArray) where {S,C,D} = CoordinateVector{S,C}(SVector{D}(v))
CoordinateVector{S,C,D}(v::AbstractArray) where {S,C,D} = CoordinateVector{S,C}(SVector{D}(v))

CoordinateVector{S,C,D,T}(t::Tuple) where {S,C,D,T} = CoordinateVector{S,C}(SVector{D,T}(t))
CoordinateVector{S,C,D,T}(v::StaticArray) where {S,C,D,T} = CoordinateVector{S,C}(SVector{D,T}(v))
CoordinateVector{S,C,D,T}(v::AbstractArray) where {S,C,D,T} = CoordinateVector{S,C}(SVector{D,T}(v))

# Indexing
Base.getindex(c::CoordinateVector, i::Int) = c.vector[i]

Base.Tuple(c::CoordinateVector) = Tuple(c.vector)

function Base.show(io::IO, c::CoordinateVector{S,C}) where {S,C}
    print(io, CoordinateVector{S,C}, '(')
    join(io, c, ", ")
    print(io, ')')
end

# Math operations
function Base.:+(c1::CoordinateVector, c2::CoordinateVector)
    require_same_space(c1, c2)
    require_same_coordinate(c1, c2)
    return promote_type(typeof(c1), typeof(c2))(c1.vector + c2.vector)
end

Base.:-(c::CoordinateVector{S,C}) where {S,C} = CoordinateVector{S,C}(-c.vector)
Base.:-(c1::CoordinateVector, c2::CoordinateVector) = +(c1, -c2)

Base.:*(s::Real, c::CoordinateVector{S,C}) where {S,C} = CoordinateVector{S,C}(s * c.vector)
Base.:*(c::CoordinateVector{S,C}, s::Real) where {S,C} = CoordinateVector{S,C}(c.vector * s)

# Dot product as multiplication
function Base.:*(c1::CoordinateVector, c2::CoordinateVector)
    require_dual_space(c1, c2)
    require_same_coordinate(c1, c2)
    return dot(c1.vector, c2.vector)
end

Base.:/(c::CoordinateVector{S,C}, s::Real) where {S,C} = CoordinateVector{S,C}(c.vector / s)
Base.:(//)(c::CoordinateVector{S,C}, s::Real) where {S,C} = CoordinateVector{S,C}(c.vector // s)

# Defined separately because s / c.vector is not a `StaticArray`
function Base.:/(s::Real, c::CoordinateVector{S,C,D}) where {S,C,D}
    return CoordinateVector{inverse_space(S),C,D}(s / c.vector)
end

Base.:\(s::Real, c::CoordinateVector{S,C}) where {S,C} = CoordinateVector{S,C}(s \ c.vector)
# Base.:\(c::CoordinateVector{S,C}, s::Real) where {S,C} = CoordinateVector{S,C}(c.vector \ s)

#---Shift vectors----------------------------------------------------------------------------------#
"""
    ShiftVector{S<:BySpace,D,T} <: StaticVector{D,T}

A vector in fractional coordinates representing a shift of a lattice or lattice dataset from the
origin. This wraps a `CoordinateVector{S,ByFractionalCoordinate,D,T}` with an optional weight
parameter that may be useful when working with symmetrical structures or k-points in the
irreducible Brillouin zone. If it is not explicitly set, it defaults to 1.

If constructors do not explicitly reference an element type, the element type is automatically
inferred by promoting the types of the arguments.

# Type aliases

`ShiftVector` is a general way of working with vectors which shift a lattice or data within it, and
for this reason we define an alias for representing k-points:

    const KPoint = ShiftVector{ByReciprocalSpace}
"""
struct ShiftVector{S<:BySpace,D,T} <: StaticVector{D,T}
    vector::CoordinateVector{S,ByFractionalCoordinate,D,T}
    weight::T
    function ShiftVector(v::CoordinateVector{S,ByFractionalCoordinate}, wt::Real = true) where S
        return new{S, length(v), promote_type(eltype(v), typeof(wt))}(v, wt)
    end
end

BySpace(::Type{<:ShiftVector{S}}) where S = S()
ByCoordinate(::Type{<:ShiftVector}) = ByFractionalCoordinate()

"""
    KPoint{D,T} (alias for ShiftVector{ByReciprocalSpace,D,T})

Represents a k-point, or an offset from the Γ point (origin) of reciprocal space, often used in band
structures, wavefunctions, or other electronic data.

For more information about this type, see [`ShiftVector`](@ref).
"""
const KPoint = ShiftVector{ByReciprocalSpace}

function ShiftVector{S}(v::StaticVector, wt::Real = true) where S 
    return ShiftVector(CoordinateVector{S,ByFractionalCoordinate}(v), wt)
end

ShiftVector{S}(::StaticArray, wt::Real = true) where S = array_not_flattened()
ShiftVector{S}(t::Tuple, wt::Real = true) where S = ShiftVector{S}(SVector(t), wt)
ShiftVector{S}(x::Real...; weight::Real = true) where S = ShiftVector{S}(SVector(x), weight)

ShiftVector{S,D}(t::Tuple, wt::Real = true) where {S,D} = ShiftVector{S}(SVector{D}(t), wt)
ShiftVector{S,D}(v::StaticArray, wt::Real = true) where {S,D} = ShiftVector{S}(SVector{D}(v), wt)
ShiftVector{S,D}(v::AbstractArray, wt::Real = true) where {S,D} = ShiftVector{S}(SVector{D}(v), wt)

ShiftVector{S,D,T}(t::Tuple, wt::Real = true) where {S,D,T} = ShiftVector{S}(SVector{D,T}(t), T(wt))

function ShiftVector{S,D,T}(v::StaticArray, wt::Real = true) where {S,D,T}
    return ShiftVector{S}(SVector{D,T}(v), T(wt))
end

function ShiftVector{S,D,T}(v::AbstractArray, wt::Real = true) where {S,D,T}
    return ShiftVector{S}(SVector{D,T}(v), T(wt))
end

# Hashing and equality
Base.hash(s::ShiftVector, h::UInt) = hash(s.vector, hash(s.weight, h))

function Base.:(==)(u::ShiftVector{S1}, v::ShiftVector{S2}) where {S1,S2}
    return (S1 === S2 && u.vector == v.vector && u.weight == v.weight)
end

# Indexing
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

Base.summary(io::IO, s::ShiftVector) = print(io, typeof(s), " with weight ", s.weight)
Base.show(io::IO, s::ShiftVector) = print(io, typeof(s), '(', s.vector, ", ", s.weight, ')')

function Base.show(io::IO, k::KPoint)
    print(io, KPoint, '(')
    join(io, k.vector, ", ")
    print(io, ", weight = ", weight(k), ')')
end
