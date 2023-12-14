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
        basis::AbstractMatrix,
        shift::AbstractVector = zero(ShiftVector{S,D,T}),
        transform::AbstractMatrix = LinearAlgebra.I
    ) where {S,D,T}
        _transform = SVector{D,T}.(eachcol(SMatrix{D,D}(transform)))
        return new(basis, shift, _transform)
    end
end

# Custom property for transform matrix
function Base.getproperty(l::LatticeDataMap, s::Symbol)
    s === :transform && return hcat(l.transform...)
    return getfield(l, s)
end

function Base.propertynames(l::LatticeDataMap; private = false)
    names = tuple(fieldnames(typeof(l))..., :transform) 
    return private ? names : names[1,2,4]
end

# Constructors with varying numbers of type parameters and arguments
function LatticeDataMap{S,D}(
    basis::AbstractMatrix,
    shift::AbstractVector = zero(ShiftVector{S,D}),
    transform::AbstractMatrix = LinearAlgebra.I
) where {S,D}
    return LatticeDataMap{S,D,promote_eltype(basis, shift)}(basis, shift, transform)
end

function LatticeDataMap{S}(
    basis::StaticMatrix{D,D},
    shift::AbstractVector = zero(ShiftVector{S,D}),
    transform::AbstractMatrix = LinearAlgebra.I
) where {S,D}
    return LatticeDataMap{S,D}(basis, shift, transform)
end

function LatticeDataMap(
    basis::LatticeBasis{S,D},
    shift::Union{CoordinateVector,ShiftVector} = zero(ShiftVector{S,D}),
    transform::AbstractMatrix = LinearAlgebra.I
) where {S,D}
    BySpace(shift) isa S || error("Cannot determine whether real or reciprocal space was intended.")
    return LatticeDataMap{S,D}(basis, shift, transform)
end
