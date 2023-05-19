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
    ByRealSpace{D}

Trait for real space data in `D` dimensions.
"""
struct ByRealSpace{D} <: DataSpace{D}
end

"""
    ByReciprocalSpace{D}

Trait for reciprocal space data in `D` dimensions.
"""
struct ByReciprocalSpace{D} <: DataSpace{D}
end

"""
    ByAtom{D}

Trait for data associated with atomic positions in a crystal.
"""
struct ByAtom{D} <: DataSpace{D}
end

"""
    Electrum.DataSpace(x) -> CrystalDataTrait

Returns a trait that determines whether a data set associated with a crystal is defined in real
space (`RealSpaceData{D}()`), reciprocal space (`ReciprocalSpaceData{D}()`), or by atomic positions
(`AtomPositionData{D}`), where `D` is the number of dimensions.
"""
DataSpace(x) = DataSpace(typeof(x))
