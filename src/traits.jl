#---Trait supertype--------------------------------------------------------------------------------#
"""
    CrystalDataTrait

Subtypes of this type are traits that may be used for dispatch.
"""
abstract type CrystalDataTrait
end

#---Association of data with spatial coordinates---------------------------------------------------#
"""
    DataSpace{D} <: CrystalDataTrait

Describes the space in which a dataset is defined, which can be real space, reciprocal space, or
data associated with individual atoms in a structure.
"""
abstract type DataSpace{D} <: CrystalDataTrait
end

"""
    ByAtom{D}

Trait for data associated with atomic positions in a crystal.
"""
struct ByAtom{D} <: DataSpace{D}
end

"""
    BySpace{D}

Supertype for the `ByRealSpace{D}` and `ByReciprocalSpace{D}` traits.
"""
abstract type BySpace{D} <: DataSpace{D}
end

"""
    ByRealSpace{D}

Trait for real space data in `D` dimensions.
"""
struct ByRealSpace{D} <: BySpace{D}
end

"""
    ByReciprocalSpace{D}

Trait for reciprocal space data in `D` dimensions.
"""
struct ByReciprocalSpace{D} <: BySpace{D}
end

"""
    Electrum.inverse_space(::Type{<:BySpace}) -> Type{<:BySpace}
    Electrum.inverse_space(::BySpace) -> BySpace

Returns the space trait dual to the given space trait. Either a type or an instance may be given.

# Examples
```julia-repl
julia> Electrum.inverse_space(Electrum.ByReciprocalSpace)
Electrum.ByRealSpace

julia> Electrum.inverse_space(Electrum.ByRealSpace{3})
Electrum.ByReciprocalSpace{3}

julia> Electrum.inverse_space(Electrum.ByRealSpace{3}())
Electrum.ByReciprocalSpace{3}()
```
"""
inverse_space(::Type{ByRealSpace}) = ByReciprocalSpace
inverse_space(::Type{ByReciprocalSpace}) = ByRealSpace
inverse_space(::Type{ByRealSpace{D}}) where D = ByReciprocalSpace{D}
inverse_space(::Type{ByReciprocalSpace{D}}) where D = ByRealSpace{D}
inverse_space(::T) where T<:BySpace = inverse_space(T)()

"""
    Electrum.DataSpace(x) -> DataSpace

Returns a trait that determines whether a data set associated with a crystal is defined in real
space (`ByRealSpace{D}()`), reciprocal space (`ByReciprocalSpace{D}()`), or by atomic positions
(`ByAtom{D}`), where `D` is the number of dimensions.

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

"""
    Electrum.dimension(::DataSpace{D}) = D
    Electrum.dimension(::Type{<:DataSpace{D}}) = D

Returns the number of dimensions (or other object representing the dimensionality) associated with a
`DataSpace` trait.
"""
dimension(::DataSpace{D}) where D = D
dimension(::Type{<:DataSpace{D}}) where D = D

"""
    Electrum.dimension(x)

Infers the number of dimensions `Electrum.DataSpace(x)` from the result of `DataSpace(x)`.
"""
dimension(x) = dimension(DataSpace(x))

#---Coordinate type--------------------------------------------------------------------------------#
"""
    ByCoordinate{D} <: CrystalDataTrait

Describes the coordinate system associated with data. This includes `ByCartesianCoordinate{D}` and
`ByFractionalCoordinate{D}`.
"""
abstract type ByCoordinate{D} <: CrystalDataTrait
end

"""
    ByCartesianCoordinate{D} <: ByCoordinate{D}

Trait type for coordinates in `D` dimensions represented in terms of an implicit orthonormal basis
in units of bohr or rad*bohr⁻¹.
"""
struct ByCartesianCoordinate{D} <: ByCoordinate{D}
end

"""
    ByFractionalCoordinate{D} <: ByCoordinate{D}

Trait type for coordinates in `D` dimensions whose values depend on a choice of basis, usually the
basis vectors describing the lattice in which the coordinate is contained.
"""
struct ByFractionalCoordinate{D} <: ByCoordinate{D}
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
