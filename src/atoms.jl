"""
    AtomPosition{D} <: AbstractRealSpaceData{D}

Contains information about an atom in a crystal, including atomic number, the name of an atom,
and its position. Atomic names may be set arbitrarily, but will default to the atomic symbol if not
provided explicitly.

An atomic number of 0 indicates a dummy atom - this does not represent a real atom, but a position 
of note within the structure. If no name is given for a dummy atom, it will be given the empty 
string.
"""
struct AtomPosition{D} <: AbstractRealSpaceData{D}
    name::String
    num::Int
    pos::SVector{D,Float64}
    function AtomPosition{D}(
        name::AbstractString,
        num::Integer,
        pos::AbstractVector{<:Real}
    ) where D
        return new{D}(name, num, pos)
    end
end

function AtomPosition{D}(num::Integer, pos::AbstractVector{<:Real}) where D
    return AtomPosition{D}(ELEMENTS[num], num, pos)
end

function AtomPosition{D}(name::AbstractString, pos::AbstractVector{<:Real}) where D
    return AtomPosition{D}(name, get(ELEMENT_LOOKUP, name, 0), pos)
end

"""
    atomicno(a::AtomPosition) -> Int

Gets the atomic number of an atom in an atomic position.
"""
atomicno(a::AtomPosition) = a.num

"""
    coord(a::AtomPosition{D}) -> SVector{D,Float64}
"""
coord(a::AtomPosition) = a.pos


"""
    AtomList{D} <: AbstractRealSpaceData{D}

A list of atomic positions with an associated basis. Atomic coordinates are given in terms of the 
specified basis. If the basis is a zero matrix, the coordinates are assumed to be given in 
angstroms.
"""
struct AtomList{D} <: AbstractRealSpaceData{D}
    basis::SMatrix{D,D,Float64}
    coord::Vector{AtomPosition{D}}
    function AtomList{D}(
        basis::AbstractMatrix{<:Number},
        coord::AbstractVector{AtomPosition{D}}
    ) where D
        # Need to check if the basis is valid (only if it's nonzero, skip otherwise)
        iszero(basis) || lattice_sanity_check(basis)
        # Remove any duplicate atoms if they come up
        return new{D}(basis, unique(coord))
    end
end

function AtomList{D}(
    basis::AbstractLattice{D},
    coord::AbstractVector{AtomPosition{D}};
    prim=false
) where D
    # Use either the primitive or conventional vectors depending on user input
    if prim
        return AtomList{D}(basis.prim, unique(coord))
    else
        return AtomList{D}(basis.conv, unique(coord))
    end
end

# Do it without adding a lattice
function AtomList{D}(coord::AbstractVector{AtomPosition{D}}) where D
    return AtomList{D}(zeros(SMatrix{D,D,Float64}), coord)
end

# Check if an AtomList{D} is empty
Base.isempty(l::AtomList{D}) where D = isempty(l.coord)

# Iterate through an AtomList
Base.iterate(l::AtomList{D}) where D = iterate(l.coord)
Base.iterate(l::AtomList{D}, n) where D = iterate(l.coord, n)

natom(l::AtomList{D}) where D = length(l.coord)

"""
    cartesian(l::AtomList{D}) -> AtomList{D}

Converts an `AtomList` from reduced coordinates (relative to some crystal basis) to Cartesian
coordinates in space.
"""
function cartesian(l::AtomList{D}) where D
    # If there are no basis vectors specified, just return the original list
    # If the basis vectors are zero we're assuming Cartesian coordinates
    l.basis == zeros(SMatrix{D,D,Float64}) && return l
    newlist = map(a -> cartesian(l.basis, a), l.coord)
    return AtomList{D}(zeros(SMatrix{D,D,Float64}), newlist)
end

"""
    cartesian(b::AbstractMatrix{<:Real}, a::AtomPosition{D}) -> AtomPosition{D}

Converts an `AtomPosition` defined in terms of basis `b` to a new `AtomPosition` defined in 
Cartesian coordinates.
"""
function cartesian(b::AbstractMatrix{<:Real}, a::AtomPosition{D}) where D
    newpos = b * a.pos
    return AtomPosition{D}(a.name, a.num, newpos)
end