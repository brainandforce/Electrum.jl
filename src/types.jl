#=
IMPORTANT NOTE:

Unfortuantely, Julia does not include any structure within abstract types. This can pose a problem,
because it's not possible to immediately determine what requirements abstract types must meet.

Therefore, any declared abstract type must contain a comment describing those requirements, either 
in the docstring or in a comment within the abstract type declaration. Of course, they are not 
strictly binding, but make sure that you know why you're breaking the rules if you choose to do so.
=#

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
    AbstractCrystalData{D}

Supertype for all datasets that may be associated with a crystal.
"""
abstract type AbstractCrystalData{D}
end

"""
    AbstractRealSpaceData{D}

Supertype for crystal data that is defined in real space. This includes atomic coordinates.

Any basis that is stored with a subtype of this type will be a `RealBasis`, and calling `basis()`
on that subtype will always return a `RealBasis.`
"""
abstract type AbstractRealSpaceData{D} <: AbstractCrystalData{D}
end

"""
    AbstractReciprocalSpaceData{D}

Supertype for crystal data that is defined in reciprocal space.

Any basis that is stored with a subtype of this type will be a `ReciprocalBasis`, and calling
`basis()` on that subtype will always return a `ReciprocalBasis.`
"""
abstract type AbstractReciprocalSpaceData{D} <: AbstractCrystalData{D}
end

"""
    AbstractHKL{D}

Supertype for crystal data stored by HKL index.
"""
abstract type AbstractHKL{D,T} <: AbstractReciprocalSpaceData{D}
end

"""
    AbstractDensityOfStates

Supertype for all density of states data.
"""
abstract type AbstractDensityOfStates
end

"""
    AbstractKPoints{D}

Supertype for sets of k-points.
"""
abstract type AbstractKPoints{D} <: AbstractReciprocalSpaceData{D}
end
