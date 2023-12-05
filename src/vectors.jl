#---Coordinate vectors-----------------------------------------------------------------------------#
"""
    CoordinateVector{S<:BySpace,C<:ByCoordinate,D,T<:Real} <: StaticVector{D,T}

Represents a spatial coordinate with space given by trait `S` and coordinate trait `C`. Coordinate
values must be subtypes of `Real`.
"""
struct CoordinateVector{S<:BySpace,C<:ByCoordinate,D,T<:Real} <: StaticVector{D,T}
    data::NTuple{D,T}
    CoordinateVector{S,C,D,T}(t::Tuple) where {S,C,D,T} = new(t)
end

Base.getindex(c::CoordinateVector, i::Int) = c.data[i]
Base.Tuple(c::CoordinateVector) = c.data

BySpace(::Type{<:CoordinateVector{S}}) where S = S()
ByCoordinate(::Type{<:CoordinateVector{<:BySpace,C}}) where C = C()

function StaticArrays.similar_type(::Type{V}, ::Type{T}, ::Size{S}) where {V<:CoordinateVector,T,S}
    if isone(length(S)) 
        return CoordinateVector{typeof(BySpace(V)), typeof(ByCoordinate(V)), only(S), T}
    else
        return SArray{Tuple{S...},T}
    end
end

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

CoordinateVector{S,C,D}(t::Tuple) where {S,C,D} = CoordinateVector{S,C,D,promote_typeof(t...)}(t)
CoordinateVector{S,C}(t::Tuple) where {S,C} = CoordinateVector{S,C,length(t)}(t)

(V::Type{<:CoordinateVector{S,C}})(x::Real...) where {S,C} = V(x)

(V::Type{<:CoordinateVector{S,C,D}})(s::StaticArray) where {S,C,D} = V(Tuple(s))
(V::Type{<:CoordinateVector{S,C,D}})(a::AbstractArray) where {S,C,D} = V(NTuple{D}(a))

CoordinateVector{S,C}(s::StaticVector) where {S,C} = CoordinateVector{S,C}(Tuple(s))
CoordinateVector{S,C}(::StaticArray) where {S,C} = array_not_flattened()

function (V::Type{<:CoordinateVector{S,C,D}})(v::CoordinateVector) where {S,C,D}
    BySpace(v) isa S || error("Cannot convert between real and reciprocal space coordinates.")
    ByCoordinate(v) isa C || error(
        "Use matrix operations with lattice basis vectors to perform this conversion:\n" *
        if C <: ByCartesianCoordinate
            "Multiply fractional coordinates by the basis to obtain Cartesian coordinates."
        elseif C <: ByFractionalCoordinate
            "Right divide Cartesian coordinates by the basis to obtain fractional coordinates."
        else
            "(No help is available for this type of transformation)"
        end
    )
    return V(Tuple(v))
end

CoordinateVector{S,C}(v::CoordinateVector) where {S,C} = CoordinateVector{S,C,length(v)}(v)
CoordinateVector{S}(v::CoordinateVector) where {S} = CoordinateVector{S,typeof(ByCoordinate(v))}(v)
CoordinateVector(v::CoordinateVector) = v

Base.convert(::Type{T}, v::CoordinateVector) where  T<:CoordinateVector = T(v)

# Pretty printing
function Base.show(io::IO, c::CoordinateVector{S,C}) where {S,C}
    print(io, typeof(c), '(')
    join(io, c, ", ")
    print(io, ')')
end

# Math operations
function Base.:+(c1::CoordinateVector, c2::CoordinateVector)
    require_same_space(c1, c2)
    require_same_coordinate(c1, c2)
    return promote_type(typeof(c1), typeof(c2))(c1.data .+ c2.data)
end

Base.:-(c::CoordinateVector{S,C}) where {S,C} = CoordinateVector{S,C}(.-c.data)
# Multiplicative operations do not guarantee the coordinate system or space are known.
Base.:*(s::Real, c::CoordinateVector) = SVector(s .* c.data)
Base.:*(c::CoordinateVector, s::Real) = SVector(c.data .* s)
Base.:/(c::CoordinateVector, s::Real) = SVector(c.data ./ s)
Base.:/(s::Real, c::CoordinateVector) = SVector(s ./ c.data)
Base.:\(s::Real, c::CoordinateVector) = SVector(s .\ c.data)
Base.:\(c::CoordinateVector, s::Real) = SVector(c.data .\ s)

# Dot product as multiplication
function Base.:*(c1::CoordinateVector, c2::CoordinateVector)
    require_dual_space(c1, c2)
    require_same_coordinate(c1, c2)
    return dot(SVector(c1.data), SVector(c2.data))
end

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
    function ShiftVector(v::CoordinateVector{S,C}, wt::Real = true) where {S,C}
        C === ByFractionalCoordinate || error("Only fractional coordinates are allowed.")
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

(T::Type{<:ShiftVector{<:BySpace}})(x::Real...; weight::Real = true) = T(SVector(x), weight)

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

function Base.show(io::IO, s::ShiftVector)
    print(io, typeof(s), '(')
    join(io, s.vector, ", ")
    print(io, ", weight = ", weight(s), ')')
end

#---Shared conversion semantics--------------------------------------------------------------------#

(T::Type{<:CoordinateVector})(v::ShiftVector) = T(v.vector)
Base.convert(::Type{T}, v::ShiftVector) where T<:CoordinateVector = T(v)
Base.convert(::Type{T}, v::CoordinateVector) where T<:ShiftVector = T(v)

#---Mapping array indices to grids-----------------------------------------------------------------#
"""
    LatticeDataMap{S,D,T}

Describes how the data of a `D`-dimensional array maps to regularly spaced positions in a
`D`-dimensional periodic lattice.

This data structure does not reference the size of any array, allowing it to be reused as-is if a
grid with different dimensions is specified along with this map.

# Examples

The most straightforward manner to do this for an array `a` is to map each `CartesianIndex` `i` to a
fractional coordinate by associating the first index (`firstindex(CartesianIndices(a))`) with the
origin, then spacing out the grid points evenly along the the directions given by each lattice basis
vector. This can be accomplished by simply calling the constructor with a basis:
```julia-repl
b = RealBasis{3}([1 0 0; 0 2 0; 0 0 3]); # face-centered orthorhombic

julia> LatticeDataMap(b)

```
In some cases (particularly in reciprocal space), you may want to use a `ShiftVector` (specifically
a `KPoint`, which aliases `ShiftVector{ByReciprocalSpace}`) to offset the array data:
```julia-repl
julia> LatticeDataMap(dual(b), KPoint(1/4, 1/4, 1/4))

```
However, the most complex case involves altering the mapping of points with a transform matrix. This
may be convenient for conventional cells whose data is provided with respect to a primitive cell. In
this case, a rational matrix transform can be used. The transform `T` defined below converts the
lattice basis vectors of a face-centered conventional cell to its primitive cell.
```julia-repl
julia> T = 1//2 * [0 1 1; 1 0 1; 1 1 0]; # gets primitive cell of face-centered cell

julia> LatticeDataMap(b, zero(ShiftVector{ByRealSpace,3}), T)

```
"""
struct LatticeDataMap{S<:BySpace,D,T}
    basis::LatticeBasis{S,D,T}
    shift::ShiftVector{S,D,T}
    _transform::SVector{D,SVector{D,Rational{Int}}}
    function LatticeDataMap{S,D,T}(
        basis::LatticeBasis,
        shift::ShiftVector = zero(ShiftVector{S,D,T}),
        transform::AbstractMatrix = LinearAlgebra.I
    ) where {S,D,T}
        _transform = SVector{D,T}.(eachcol(SMatrix{D,D}(transform)))
        return new(basis, shift, _transform)
    end
end

# Constructors with varying numbers of type parameters and arguments
function LatticeDataMap{S,D}(
    basis::LatticeBasis,
    shift::ShiftVector = zero(ShiftVector{S,D}),
    transform::AbstractMatrix = LinearAlgebra.I
) where {S,D}
    T = promote_type(eltype(basis), eltype(shift))
    return LatticeDataMap{S,D,T}(basis, shift, transform)
end

function LatticeDataMap{S}(
    basis::LatticeBasis{<:BySpace,D},
    shift::ShiftVector{<:BySpace,D} = zero(ShiftVector{S,D}),
    transform::AbstractMatrix = LinearAlgebra.I
) where {S,D}
    return LatticeDataMap{S,D}(basis, shift, transform)
end

function LatticeDataMap(
    basis::LatticeBasis{S1,D},
    shift::ShiftVector{S2,D} = zero(ShiftVector{S1,D}),
    transform::AbstractMatrix = LinearAlgebra.I
) where {S1,S2,D}
    S1 === S2 || error("Cannot determine whether real or reciprocal space was intended.")
    return LatticeDataMap{S1,D}(basis, shift, transform)
end

#---Data defined on regularly spaced grids on lattices---------------------------------------------#
"""
    LatticeData{S,D,M<:LatticeDataMap{C,D},T,A<:AbstractArray{T,D}} <: AbstractArray{T,D}

Maps an array of type `A` onto regularly spaced points of a lattice with a mapping of type `M`.

# Indexing

Arrays wrapped by this data structure become zero-based arrays with periodic Cartesian indexing. For
this reason, all indexing operations are guaranteed to be inbounds.
"""
struct LatticeData{S,D,M<:LatticeDataMap{S,D},T,A<:AbstractArray{T,D}} <: AbstractArray{T,D}
    data::A
    map::M
end

Base.size(l::LatticeData) = size(l.data)
Base.axes(l::LatticeData) = map(x -> x .- first.(axes(l.data)), axes(l.data))

# Index style matches that of backing array
Base.IndexStyle(::Type{<:LatticeData{<:Any,<:Any,<:Any,<:Any,A}}) where A = IndexStyle(A)

"""
    Electrum._to_array_index(l::LatticeData, dim::Integer, i)

Converts an index `i` into a valid index of the backing array of `l` along dimension `dim`. To
convert an arbitrary integer index to a valid index of the backing array, it uses the formula

    mod(i, size(l, dim)) + first(axes(l.data, dim))

For `CartesianIndex` arguments, `dims` refers to the first dimension indexed by `i`, and the
boundary conditions for the tail indices are automatically calculated.

This formula is extended to `AbstractArray` arguments, and colon arguments are unchanged.
"""
function _to_array_index(l::LatticeData, dim::Integer, i::Integer)
    return Int(mod(i, size(l, dim)) + first(axes(l.data, dim)))
end

function _to_array_index(l::LatticeData, dim::Integer, i::CartesianIndex)
    dimrange = (dim - 1) .+ 1:length(i)
    return CartesianIndex(mod.(Tuple(i), size(l)[dimrange]) .+ first.(axes(l.data))[dimrange])
end

function _to_array_index(l::LatticeData, dim::Integer, i::AbstractArray{<:Integer})
    return Int.(mod.(i, size(l, dim)) .+ first(axes(l.data, dim)))
end

function _to_array_index(l::LatticeData, dim::Integer, i::CartesianIndices)
    dimrange = (dim - 1) .+ 1:length(i.indices)
    return CartesianIndices(mod.(i.indices, size(l)[dimrange] .+ first.(axes(l.data)[dimrange])))
end

_to_array_index(::LatticeData, ::Integer, ::Colon) = Colon()

"""
    Electrum._meta_indices(I::Tuple)

Returns a `NTuple{length(I),Int}` corresponding to the dimensions indexed by `I`. This is used to
correct for the presence of any `CartesianIndex` or `CartesianIndices` arguments which may index
multiple array dimensions.
"""
function _meta_indices(I::Tuple)
    v = Tuple(eachindex(I))
    extra_length = 0
    for (n, x) in enumerate(I)
        v = Base.setindex(v, v[n] + extra_length, n)
        if x isa CartesianIndex
            extra_length += length(x) - 1
        elseif x isa AbstractArray{CartesianIndex{D},D} where D
            extra_length += ndims(x) - 1
        end
    end
    return v
end

"""
    Electrum._to_array_indices(l::LatticeData, I::Tuple)

Converts a set of indices from a `getindex` or `setindex!` call to valid indices of the backing
array of `l`.
"""
function _to_array_indices(l::LatticeData, I::Tuple)
    return map((x, dim) -> _to_array_index(l, dim, x), I, _meta_indices(I))
end

Base.getindex(l::LatticeData, i...) = @inbounds getindex(l.data, _to_array_indices(l, i)...)
Base.setindex!(l::LatticeData, x, i...) = @inbounds setindex!(l.data, x, _to_array_indices(l, i)...)
