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
    function AtomPosition(
        name::AbstractString,
        num::Integer,
        pos::SVector{D,<:Real}
    ) where D
        # Don't allow empty names when an atomic number is passed.
        name = isempty(name) && num > 0 ? ELEMENTS[num] : name
        return new{D}(name, num, pos)
    end
end

function AtomPosition(name::AbstractString, pos::SVector{D,<:Real}) where D
    return AtomPosition(name, get(ELEMENT_LOOKUP, name, 0), pos)
end

function AtomPosition(num::Integer, pos::SVector{D,<:Real}) where D
    return AtomPosition(ELEMENTS[num], num, pos)
end

# The methods below require a type parameter in the constructor

function AtomPosition{D}(name::AbstractString, num::Integer, pos::AbstractVector{<:Real}) where D
    return AtomPosition(name, num, SVector{D,Float64}(pos))
end

function AtomPosition{D}(name::AbstractString, pos::AbstractVector{<:Real}) where D
    return AtomPosition(name, SVector{D,Float64}(pos))
end

function AtomPosition{D}(num::Integer, pos::AbstractVector{<:Real}) where D
    return AtomPosition(num, SVector{D,Float64}(pos))
end

"""
    atomname(a::AtomPosition) -> String

Returns the name of the atom. Note that this returns the name that was provided during
construction, which is not guaranteed to be the same as the name of the element that has the given
atomic number.
"""
atomname(a::AtomPosition) = a.name

"""
    atomicno(a::AtomPosition) -> Int

Gets the atomic number of an atom in an atomic position.
"""
atomicno(a::AtomPosition) = a.num

"""
    coord(a::AtomPosition{D}) -> SVector{D,Float64}

Returns the coordinate associated with an atomic position.
"""
coord(a::AtomPosition) = a.pos

Base.getindex(a::AtomPosition, ind) = a.pos[ind]

"""
    AtomList{D} <: AbstractRealSpaceData{D}

A list of atomic positions with an associated basis. Atomic coordinates are given in terms of the 
specified basis. If the basis is a zero matrix, the coordinates are assumed to be given in 
angstroms.
"""
struct AtomList{D} <: AbstractRealSpaceData{D}
    basis::BasisVectors{D}
    coord::Vector{AtomPosition{D}}
    function AtomList(
        basis::BasisVectors{D},
        coord::AbstractVector{<:AtomPosition{D}}
    ) where D
        # Remove any duplicate atoms if they come up
        return new{D}(basis, unique(coord))
    end
end

function AtomList(
    basis::AbstractLattice{D},
    coord::AbstractVector{<:AtomPosition{D}};
    prim=false
) where D
    # Use either the primitive or conventional vectors depending on user input
    if prim
        return AtomList(basis.prim, unique(coord))
    else
        return AtomList(basis.conv, unique(coord))
    end
end

# Do it without adding a lattice
function AtomList(coord::AbstractVector{AtomPosition{D}}) where D
    return AtomList(zeros(BasisVectors{D}), coord)
end

# Check if an AtomList{D} is empty
Base.isempty(l::AtomList) = isempty(l.coord)
Base.getindex(l::AtomList, ind) = l.coord[ind]
Base.length(l::AtomList) = length(l.coord)
Base.size(l::AtomList) = size(l.coord)

# Iterate through an AtomList
Base.iterate(l::AtomList) = iterate(l.coord)
Base.iterate(l::AtomList, n) = iterate(l.coord, n)

basis(l::AtomList) = l.basis
coord(l::AtomList) = l.coord
natom(l::AtomList) = length(coord(l))

"""
    cartesian(l::AtomList{D}) -> AtomList{D}

Converts an `AtomList` from reduced coordinates (relative to some crystal basis) to Cartesian
coordinates in space.
"""
function cartesian(l::AtomList{D}) where D
    # If there are no basis vectors specified, just return the original list
    # If the basis vectors are zero we're assuming Cartesian coordinates
    basis(l) == zeros(BasisVectors{D}) && return l
    newlist = map(a -> cartesian(l.basis, a), l.coord)
    return AtomList(zeros(BasisVectors{D}), newlist)
end

# TODO: clean up these method definitions into something that makes more sense
"""
    cartesian(b::AbstractMatrix{<:Real}, a::AtomPosition{D}) -> AtomPosition{D}

Converts an `AtomPosition` defined in terms of basis `b` to a new `AtomPosition` defined in 
Cartesian coordinates.
"""
function cartesian(b::AbstractMatrix{<:Real}, a::AtomPosition{D}) where D
    newpos = b * a.pos
    return AtomPosition{D}(a.name, a.num, newpos)
end

cartesian(b::BasisVectors{D}, a::AtomPosition{D}) where D = cartesian(matrix(b), a)

"""
    reduce_coords(basis::AbstractMatrix{<:Real}, a::AtomPosition; incell=false)

Convert a coordinate from a Cartesian basis to the crystal basis. If `incell` is true, the position
vector components will be truncated so they lie within the cell bounds (between 0 and 1).
"""
function reduce_coords(
    basis::AbstractMatrix{<:Real},
    a::AtomPosition{D};
    incell::Bool=false
) where D
    v = basis\coord(a)
    # If it needs to be in a cell, truncate the integer portion of the elements 
    v = v - incell * floor.(v)
    return AtomPosition{D}(name(a), atomicno(a), v)
end

"""
    reduce_coords(basis::AbstractMatrix{<:Real}, a::AbstractVector{AtomPosition}; incell=false)

Convert a vector of atomic coordinates from a Cartesian basis to the crystal basis. If `incell` is
true, the position vector components will be truncated so they lie within the cell bounds (between
0 and 1).
"""
function reduce_coords(
    basis::AbstractMatrix{<:Real},
    va::AbstractVector{AtomPosition{D}};
    incell::Bool = false
) where D
    va_new = map(va) do a
        reduce_coords(basis, a, incell = incell)
    end
    return AtomPosition{D}(basis, va_new)
end

"""
    supercell(l::AtomList, M::AbstractMatrix{<:Integer}) -> AtomList

Creates a new `AtomList` with the basis vectors of a supercell generated by transforming the basis 
vectors of the space by `M`. This will also generate new atomic positions to fill the cell.
"""
function supercell(l::AtomList{D}, M::AbstractMatrix{<:Integer}) where D
    # Convert the provided basis vectors to the supercell BasisVectors
    # Conversion to upper triangular form should make this easier to work with
    scb = triangularize(basis(l), M)
    # Use LU decomposition to generate the translation bounds
    # Convert to a positive integer to generate valid indices
    tmax = SVector{D,Int}(abs.(diag(lu(M).U)))
    # Generate all the sites in the supercell where new atoms have to be placed
    newpts = [(v.I .- 1) ./ tmax for v in CartesianIndices(Tuple(tmax))]
    # Shift everything over for each new point
    sclist = [
        AtomPosition(
            atomname(atom),
            atomicno(atom),
            (v = M \ coord(atom) + d;  v - floor.(v)) # Keep the new atoms inside the supercell
        )
        for atom in l.coord, d in newpts
    ]
    return AtomList(scb, vec(sclist))
end

supercell(l::AtomList, v::AbstractVector{<:Integer}) = supercell(l, diagm(v))

"""
    remove_dummies(l::AtomList) -> AtomList

Removes dummy atoms from an `AtomList`.
"""
 remove_dummies(l::AtomList) = AtomList(basis(l), filter(x -> !iszero(atomicno(x)), l.coord))

"""
    natomtypes(l::AtomList; dummy=false) -> Int

Returns the number of types of atoms. The `dummy` keyword controls whether dummy atoms are counted
as a separate atom type (`false` by default).
"""
function natomtypes(l::AtomList; dummy=false)
    # Collect all of the atom types
    types = [atomicno(a) for a in l]
    # Subtract any dummy atoms
    return length(unique(types)) - (!dummy * (0 in types))
    # (There's probably a more efficient implementation, but this works)
end

"""
    atomnames(l::AtomList; dummy=false)

Returns the names of all the atoms in `l`. By default, dummy atoms are not included, but this may
be changed by setting `dummy=true`.
"""
function atomnames(l::AtomList; dummy::Bool=false)
    # Store results here
    names = String[]
    for atom in l
        # Only add dummy atoms if requested
        if atomicno(atom) != 0 || dummy
            push!(names, atomname(atom))
        end
    end
    return unique(names)
end
