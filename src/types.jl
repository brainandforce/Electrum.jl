#=
IMPORTANT NOTE:

Unfortuantely, Julia does not include any structure within abstract types. This can pose a problem,
because it's not possible to immediately determine what requirements abstract types must meet.

Therefore, any declared abstract type must contain a comment describing those requirements, either 
in the docstring or in a comment within the abstract type declaration. Of course, they are not 
strictly binding, but make sure that you know why you're breaking the rules if you choose to do so.
=#

"""
    AbstractAtomList{D}

Supertype for lists of atomic positions in `D` dimensions.
"""
abstract type AbstractAtomList{D}
end

"""
    AbstractAtomPosition{D}

Supertype that describes atomic positions in `D` dimensions, which include name, coordinate, and
occupancy information.
"""
abstract type AbstractAtomPosition{D}
end

"""
    AbstractBasis{D}

Supertype for sets of basis vectors in `D` dimensions.

This supertype includes the `RealBasis{D}` and `ReciprocalBasis{D}` types, which explicitly 
indicate their units (assumed to be either angstroms or inverse angstroms).

Members of `AbstractBasis` must implement the following checks:
  * That the basis vectors are linearly independent and form a right-handed coordinate system, 
unless an explicit zero basis is constructed (implying no periodicity).
"""
abstract type AbstractBasis{D}
end

"""
    AbstractCrystal{D}

A crystal structure in `D` dimensions, containing information about the lattice, space group, and
atoms contained within the crystal.
"""
abstract type AbstractCrystal{D}
end

"""
    AbstractDataGrid{D,T}

Supertype for crystal data associated with a grid of entries of type `T` in real or reciprocal
space of dimension `D`.
"""
abstract type AbstractDataGrid{D,T} <: AbstractArray{T,D}
end

"""
    AbstractDensityOfStates

Supertype for all density of states data.
"""
abstract type AbstractDensityOfStates
end

"""
    AbstractKPointSet{D}

Supertype for sets of k-points in `D` dimensions, either provided as explicit lists or as a
generator (such as a matrix defining a mesh).
"""
abstract type AbstractKPointSet{D}
end

#---Traits----------------------------------------------------------------------------------------#

"""
    CrystalDataTraits

Subtypes of this type are traits that may be used for dispatch.
"""
abstract type CrystalDataTrait
end

"""
    RealSpaceData{D}

Trait for real space data in `D` dimensions.
"""
struct RealSpaceData{D} <: CrystalDataTrait
end

"""
    ReciprocalSpaceData{D}

Trait for reciprocal space data in `D` dimensions.
"""
struct ReciprocalSpaceData{D} <: CrystalDataTrait
end

"""
    AtomPositionData{D}

Trait for data associated with atomic positions in a crystal.
"""
struct AtomPositionData{D} <: CrystalDataTrait
end

"""
    Electrum.data_space(x) -> CrystalDataTrait

Returns a trait that determines whether a data set associated with a crystal is defined in real
space (`RealSpaceData{D}()`), reciprocal space (`ReciprocalSpaceData{D}()`), or by atomic positions
(`AtomPositionData{D}`), where `D` is the number of dimensions.
"""
function data_space end
