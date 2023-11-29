#---Traits----------------------------------------------------------------------------------------#
"""
    CrystalDataTrait

Subtypes of this type are traits that may be used for dispatch.
"""
abstract type CrystalDataTrait
end

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
