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
    basis::RealBasis{D}
    coord::Vector{AtomPosition{D}}
    function AtomList(
        basis::AbstractBasis{D},
        atoms::AbstractVector{<:AtomPosition{D}}
    ) where D
        # Warn on atomic pairs that are *really* close
        # Transforming the data to Cartesian coordinates isn't necessary
        for m in eachindex(atoms)
            for n in m+1:length(atoms)
                (a,b) = atoms[[m,n]]
                dist = norm(coord(a) - coord(b))
                sametype = atomname(a) == atomname(b) && atomicno(a) == atomicno(b)
                if isapprox(dist, 0, atol=sqrt(eps(Float64)))
                    if sametype
                        @warn string(
                            "Atoms $m and $n are very close or identical!\n",
                            "Atom $m: ", repr("text/plain", a),
                            "Atom $n: ", repr("text/plain", b),
                            "Distance (in supplied basis): $dist"
                        )
                    else
                        @warn string(
                            "Atoms $m and $n are very close! Is this a mixed occupancy site?\n",
                            "Atom $m: ", repr("text/plain", a),
                            "Atom $n: ", repr("text/plain", b),
                            "Distance (in supplied basis): $dist"
                        )
                    end
                end
            end
        end
        # Remove any duplicate atoms if they come up
        return new{D}(basis, unique(atoms))
    end
end

# Do it without adding a lattice
function AtomList(coord::AbstractVector{AtomPosition{D}}) where D
    return AtomList(zeros(RealBasis{D}), coord)
end

# Check if an AtomList{D} is empty
Base.isempty(l::AtomList) = isempty(l.coord)
Base.firstindex(l::AtomList) = firstindex(l.coord)
Base.lastindex(l::AtomList) = lastindex(l.coord)
Base.getindex(l::AtomList, ind) = l.coord[ind]
Base.length(l::AtomList) = length(l.coord)
Base.size(l::AtomList) = size(l.coord)

# Iterate through an AtomList
Base.iterate(l::AtomList) = iterate(l.coord)
Base.iterate(l::AtomList, n) = iterate(l.coord, n)

Base.filter(f, l::AtomList) = AtomList(basis(l), filter(f, l.coord))
# More filter definitons that might be useful
"""
    Base.filter(name::AbstractString, l::AtomList) -> AtomList
    Base.filter(no::Integer, l::AtomList) -> AtomList

Filters an `AtomList` by atomic name or number.
"""
Base.filter(name::AbstractString, l::AtomList) = filter(a -> atomname(a) == name, l)
Base.filter(no::Integer, l::AtomList) = filter(a -> atomicno(a) == no, l)

"""
    sort_atomicno(l::AtomList; rev=false) -> AtomList

Sorts an `AtomList` by atomic number, return a new `AtomList` with the sorted elements.

The `rev` keyword is supported if a reversed order is desired.
"""
sort_atomicno(l::AtomList; kwargs...) = AtomList(basis(l), sort(l.coord, by=atomicno; kwargs...))

basis(l::AtomList) = l.basis
coord(l::AtomList) = l.coord
"""
    natom(l::AtomList) -> Int

Gets the number of atoms in an `AtomList`.
"""
natom(l::AtomList) = length(coord(l))

"""
    cartesian(l::AtomList{D}) -> AtomList{D}

Converts an `AtomList` from reduced coordinates (relative to some crystal basis) to Cartesian
coordinates in space.
"""
function cartesian(l::AtomList{D}) where D
    # If there are no basis vectors specified, just return the original list
    # If the basis vectors are zero we're assuming Cartesian coordinates
    basis(l) == zeros(RealBasis{D}) && return l
    newlist = map(a -> cartesian(l.basis, a), l.coord)
    return AtomList(zeros(RealBasis{D}), newlist)
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

cartesian(b::AbstractBasis{D}, a::AtomPosition{D}) where D = cartesian(matrix(b), a)

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
    deduplicate(l::AbstractVector{<:AtomPosition}) -> Vector{<:AtomPosition}

Removes atoms that are duplicates or close to being duplicates. In order to be considered
duplicates, the atoms must have both the same name and atomic number, and their coordinates must
be approximately equal (to within a total distance of `sqrt(eps(Float64))`).
"""
function deduplicate(l::AbstractVector{<:AtomPosition})
    # Store the indices of the atoms to keep in a vector
    kept_inds = collect(eachindex(l))
    # Loop through each pair of atoms
    for m in eachindex(l)
        for n in m+1:length(l)
            (a,b) = l[[m,n]]
            if atomname(a) == atomname(b) && atomicno(a) == atomicno(b)
                # If they're too close, set the kept_inds entry to zero
                if isapprox(norm(coord(a)), norm(coord(b)), atol=sqrt(eps(Float64)))
                    kept_inds[n] = 0
                    @debug string(
                        "Removing atom $n:\n",
                        "Atom $m: ", repr("text/plain", a),
                        "Atom $n: ", repr("text/plain", b)
                    )
                end
            end
        end
    end
    # Filter kept_inds of zeros to get the list of kept atoms
    return l[filter(!iszero, kept_inds)]
end

"""
    deduplicate(l::AtomList) -> AtomList

Places all atoms inside the cell by truncating the integer portion of the coordinate. Any atoms
of the same type (meaning same atomic number and name) that are an extremely small distance apart
will be removed.
"""
deduplicate(l::AtomList) = iszero(basis(l)) ? l : AtomList(basis(l), deduplicate(l.coord))

"""
    supercell(l::AtomList, M::AbstractMatrix{<:Integer}) -> AtomList

Creates a new `AtomList` with the basis vectors of a supercell generated by transforming the basis 
vectors of the space by `M`. This will also generate new atomic positions to fill the cell.
"""
function supercell(l::AtomList{D}, M::AbstractMatrix{<:Integer}) where D
    # Convert the provided basis vectors to the supercell basis
    # Conversion to upper triangular form should make this easier to work with
    scb = triangularize(basis(l), M)
    # Use LU decomposition to generate the translation bounds (diagonal of U matrix)
    decomp = lu(M)
    # Convert to a positive integer to generate valid indices
    # Make sure that the matrix is correctly permuted for this process
    tmax = SVector{D,Int}(abs.(diag(decomp.U))[decomp.p])
    # Generate all the sites in the supercell where new atoms have to be placed
    newpts = vec([SVector(v.I .- 1) for v in CartesianIndices(Tuple(tmax))])
    # Move all positions into the cell; remove duplciates
    dedup = deduplicate(
        [AtomPosition(atomname(a), atomicno(a), mod.(coord(a), 1)) for a in l.coord]
    )
    # Shift everything over for each new point
    sclist = [
        AtomPosition(
            atomname(atom),
            atomicno(atom),
            # Keep the new atoms inside the supercell
            mod.(SMatrix{D,D}(M)\(coord(atom) + d), 1)
        )
        for atom in dedup, d in newpts
    ]
    @debug "Vector of AtomPositions to be added:\n" * repr(vec(sclist))
    return AtomList(scb, vec(sclist))
end

supercell(l::AtomList, v::AbstractVector{<:Integer}) = supercell(l, diagm(v))

"""
    remove_dummies(l::AtomList) -> AtomList

Removes dummy atoms from an `AtomList`.
"""
remove_dummies(l::AtomList) = AtomList(basis(l), filter(x -> !iszero(atomicno(x)), l.coord))

 """
    atomtypes(l::AtomList; dummy=false) -> Vector{Int}
    atomtypes(xtal::AbstractCrystal; dummy=false) -> Vector{Int}

Returns the atomic numbers of all the atoms in the `AtomList`. The `dummy` keyword controls whether
dummy atoms are counted as a separate atom type (`false` by default).
"""
function atomtypes(l::AtomList; dummy=false)
    m = dummy ? remove_dummies(l) : l
    return unique([atomicno(a) for a in m])
end

"""
    natomtypes(l::AtomList; dummy=false) -> Int
    natomtypes(xtal::AbstractCrystal; dummy=false) -> Int

Returns the number of types of atoms. The `dummy` keyword controls whether dummy atoms are counted
as a separate atom type (`false` by default).
"""
natomtypes(l::AtomList; dummy=false) = length(atomtypes(l, dummy=dummy))

"""
    atomnames(l::AtomList; dummy=false) -> Vector{String}

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

"""
    rotate(l::AtomList, M::AbstractMatrix{<:Real}, ctr::AbstractVector{<:Real}) -> AtomList

Rotates the atoms in an `AtomList` using the rotation matrix `M` applied at the position `ctr`
(specified in Cartesian coordinates).
"""
function rotate(l::AtomList, M::AbstractMatrix{<:Real}, ctr::AbstractVector{<:Real})
    @assert M' ≈ inv(M) "The supplied matrix is not orthogonal!"
    va = iszero(basis(l)) ? l.coord : cartesian(l).coord
    positions = [coord(a) - ctr for a in va]
    new_positions = [M*v + ctr for v in positions]
    return AtomList(
        [AtomPosition(atomname(a), atomicno(a), v) for (a,v) in zip(va, new_positions)]
    )
end

function rotate(l::AtomList, M::AbstractMatrix{<:Real})
    @assert M' ≈ inv(M) "The supplied matrix is not orthogonal!"
    va = iszero(basis(l)) ? l.coord : cartesian(l).coord
    new_positions = [M*v for v in coord.(va)]
    return AtomList(
        [AtomPosition(atomname(a), atomicno(a), v) for (a,v) in zip(va, new_positions)]
    )
end
