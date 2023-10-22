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
abstract type BySpace{D}
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
