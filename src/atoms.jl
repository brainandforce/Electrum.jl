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
    AtomList{D} <: AbstractRealSpaceData{D}

A list of atomic positions with an associated basis. Atomic coordinates are given in terms of the 
specified basis. If the basis is a zero matrix, the coordinates are assumed to be given in 
angstroms.
"""
struct AtomList{D} <: AbstractRealSpaceData{D}
    basis::SMatrix{D,D,Float64}
    coord::Vector{AtomPosition{D}}
end

function AtomList{D}(
    basis::AbstractMatrix{<:Number},
    coord::AbstractVector{AtomPosition{D}}
) where D
    # Need to check if the basis is valid (only if it's nonzero, skip otherwise)
    iszero(basis) || lattice_sanity_check(basis)
    # Remove any duplicate atoms if they come up
    return AtomList{D}(basis, unique(coord))
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