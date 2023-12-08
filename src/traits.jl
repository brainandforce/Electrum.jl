#---Trait supertype--------------------------------------------------------------------------------#
"""
    Electrum.CrystalDataTrait

Subtypes of this type are traits that may be used for dispatch.
"""
abstract type CrystalDataTrait
end

#---Association of data with spatial coordinates---------------------------------------------------#
"""
    Electrum.DataSpace <: Electrum.CrystalDataTrait

Describes the space in which a dataset is defined, which can be real space, reciprocal space, or
data associated with individual atoms in a structure.
"""
abstract type DataSpace <: CrystalDataTrait
end

"""
    Electrum.ByAtom <: Electrum.DataSpace

Trait for data associated with atoms in a crystal.
"""
struct ByAtom <: DataSpace
end

"""
    BySpace <: Electrum.DataSpace

Supertype for the `ByRealSpace` and `ByReciprocalSpace` traits.
"""
abstract type BySpace <: DataSpace
end

"""
    ByRealSpace <: BySpace

Trait for real space data.
"""
struct ByRealSpace <: BySpace
end

"""
    ByReciprocalSpace <: BySpace

Trait for reciprocal space data.
"""
struct ByReciprocalSpace <: BySpace
end

"""
    inv(::Type{<:BySpace}) -> Type{<:BySpace}
    inv(::BySpace) -> BySpace

Returns the space trait dual to the given space trait. Either a type or an instance may be given.

# Examples
```julia-repl
julia> inv(ByReciprocalSpace)
ByRealSpace

julia> inv(ByRealSpace())
ByReciprocalSpace()
```
"""
Base.inv(::Type{ByRealSpace}) = ByReciprocalSpace
Base.inv(::Type{ByReciprocalSpace}) = ByRealSpace
Base.inv(::T) where T<:BySpace = inv(T)()

#= TODO: resurrect this function at some point
"""
    Electrum.DataSpace(x) -> Electrum.DataSpace

Returns a trait that determines whether a data set associated with a crystal is defined in real
space (`ByRealSpace`), reciprocal space (`ByReciprocalSpace`), or by atomic positions (`ByAtom`).

By default, `DataSpace(x)` will infer the appropriate trait from the lattice basis vectors
included in `x`, assumed to be in a field named `basis`. The fallback definition is:

    DataSpace(x) = DataSpace(typeof(x)) # basis(x) falls back to x.basis
    DataSpace(::Type{T}) where T = DataSpace(fieldtype(T, :basis))

For types `T` which use a `DataSpace` trait, but do not contain a set of lattice basis vectors which
are `Electrum.LatticeBasis` objects stored in the `basis` field, it will be necessary to define
`DataSpace(::Type{T})` explicitly for that type.
"""
DataSpace(x) = DataSpace(typeof(x))
DataSpace(T::Type) = DataSpace(fieldtype(T, :basis))
=#

"""
    BySpace(T::Type) -> BySpace
    BySpace(x) -> BySpace

Returns a trait object that describes whether data is associated with real space
([`ByRealSpace`](@ref)) or reciprocal space ([`ByReciprocalSpace`](@ref)).

New types `T` should implement the `BySpace` trait by implementing `BySpace(::Type{T})`, as the
method for type instances is defined automatically.
"""
BySpace(::Type) = error("No BySpace trait is associated with this type.")
BySpace(x) = BySpace(typeof(x))

"""
    Electrum.require_same_space(args...; msg)

Checks that the `BySpace` traits of all arguments are identical, throwing an error with optional
message `msg` if this is not the case.
"""
function require_same_space(
    traits::Tuple{Vararg{BySpace}};
    msg = "Inputs must describe the same space (real or reciprocal)."
)
    return traits isa NTuple ? nothing : error(msg)
end

require_same_space(args...; msg...) = require_same_space(BySpace.(args); msg...)

"""
    Electrum.require_dual_space(a, b; msg)

Checks that the `BySpace` trait of `a` is the inverse of the `BySpace` trait of `b`,  throwing an
error with optional message `msg` if this is not the case.
"""
function require_dual_space(a, b; msg = "One input must describe the inverse space of the other.")
    return BySpace(a) === inv(BySpace(b)) ? nothing : error(msg)
end

#---Coordinate type--------------------------------------------------------------------------------#
"""
    ByCoordinate <: Electrum.CrystalDataTrait

Describes the coordinate system associated with data. This includes `ByCartesianCoordinate` and
`ByFractionalCoordinate`.
"""
abstract type ByCoordinate <: CrystalDataTrait
end

"""
    ByCartesianCoordinate <: ByCoordinate

Trait type for coordinates represented in terms of an implicit orthonormal basis with units of bohr
or rad*bohr⁻¹.
"""
struct ByCartesianCoordinate <: ByCoordinate
end

"""
    ByFractionalCoordinate <: ByCoordinate

Trait type for coordinates in `D` dimensions whose values depend on a choice of basis, usually the
basis vectors describing the lattice in which the coordinate is contained.
"""
struct ByFractionalCoordinate <: ByCoordinate
end

"""
    ByCoordinate(T::Type)
    ByCoordinate(x) = ByCoordinate(typeof(x))

Returns a trait instance corresponding to the coordinate system associated with type `T` or an
instance of that type.

`ByCoordinate(x)` when `x` is not a `Type` calls `ByCoordinate(typeof(x))`. To implement the trait
for a custom type `T`, define `ByCoordinate(::Type{T}) where T`.
"""
ByCoordinate(x) = ByCoordinate(typeof(x))
ByCoordinate(::Type) = error("No coordinate trait is defined for this type.")

"""
    Electrum.require_same_coordinate(args...; msg)

Checks that the `ByCoordinate` traits of all arguments are identical, throwing an error with 
optional message `msg` if this is not the case.
"""
function require_same_coordinate(
    traits::Tuple{Vararg{ByCoordinate}};
    msg = "Inputs must share a common coordinate system."
)
    return traits isa NTuple ? nothing : error(msg)
end

require_same_coordinate(args...; msg...) = require_same_coordinate(ByCoordinate.(args); msg...)
