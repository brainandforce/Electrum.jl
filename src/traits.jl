#---Trait supertype--------------------------------------------------------------------------------#
"""
    CrystalDataTrait

Subtypes of this type are traits that may be used for dispatch.
"""
abstract type CrystalDataTrait
end

#---Association of data with spatial coordinates---------------------------------------------------#
"""
    DataSpace <: CrystalDataTrait

Describes the space in which a dataset is defined, which can be real space, reciprocal space, or
data associated with individual atoms in a structure.
"""
abstract type DataSpace <: CrystalDataTrait
end

"""
    ByAtom

Trait for data associated with atomic positions in a crystal.
"""
struct ByAtom <: DataSpace
end

"""
    BySpace

Supertype for the `ByRealSpace` and `ByReciprocalSpace` traits.
"""
abstract type BySpace <: DataSpace
end

"""
    ByRealSpace

Trait for real space data.
"""
struct ByRealSpace <: BySpace
end

"""
    ByReciprocalSpace

Trait for reciprocal space data.
"""
struct ByReciprocalSpace <: BySpace
end

"""
    Electrum.inverse_space(::Type{<:BySpace}) -> Type{<:BySpace}
    Electrum.inverse_space(::BySpace) -> BySpace

Returns the space trait dual to the given space trait. Either a type or an instance may be given.

# Examples
```julia-repl
julia> Electrum.inverse_space(Electrum.ByReciprocalSpace)
Electrum.ByRealSpace

julia> Electrum.inverse_space(Electrum.ByRealSpace())
Electrum.ByReciprocalSpace()
```
"""
inverse_space(::Type{ByRealSpace}) = ByReciprocalSpace
inverse_space(::Type{ByReciprocalSpace}) = ByRealSpace
inverse_space(::T) where T<:BySpace = inverse_space(T)()

"""
    Electrum.DataSpace(x) -> DataSpace

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

#---Coordinate type--------------------------------------------------------------------------------#
"""
    ByCoordinate <: CrystalDataTrait

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
