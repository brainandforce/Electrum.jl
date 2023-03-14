"""
    NamedAtom

Stores information about an atom, which includes a name which may be up to 15 codepoints long, and
the atomic number.
"""
struct NamedAtom
    name::InlineString15
    num::Int
    NamedAtom(name::AbstractString, num::Integer) = new(name, num)
end

NamedAtom(num::Integer) = num in 1:118 ? NamedAtom(ELEMENTS[num], num) : NamedAtom("dummy", num)
NamedAtom(atomname::AbstractString) = NamedAtom(atomname, get(ELEMENT_LOOKUP, atomname, 0))

Base.show(io::IO, atom::NamedAtom) = print(io, "NamedAtom(\"", atom.name, "\", ", atom.num, ")")
Base.isless(a1::NamedAtom, a2::NamedAtom) = isless(a1.num, a2.num) || isless(a1.name, a2.name)

name(a::NamedAtom) = a.name
atomic_number(a::NamedAtom) = a.num
isdummy(a::NamedAtom) = (a.num == 0)

Base.convert(T::Type{<:AbstractString}, a::NamedAtom) = convert(T, name(a))
Base.convert(T::Type{<:Number}, a::NamedAtom) = convert(T, atomic_number(a))
Base.convert(::Type{NamedAtom}, x) = NamedAtom(x)

"""
    Electrum.reset_name(a::NamedAtom)

Creates a new `NamedAtom` with the same atomic number as `a` and the atom's normal symbol.
"""
reset_name(a::NamedAtom) = NamedAtom(atomic_number(a))

#---Atomic position data---------------------------------------------------------------------------#
"""
    FractionalAtomPosition{D}

Describes an atomic position within a crystal or other periodic structure. The coordinate in the
`pos` field is assumed to be given relative to the basis vectors of the structure.

Occupancy information is provided in the `occ` field. Note that no checking is done to ensure that
the occupancy is a reasonable value.
"""
struct FractionalAtomPosition{D} <: AbstractAtomPosition{D}
    atom::NamedAtom
    pos::SVector{D,Float64}
    occ::Float64
    function FractionalAtomPosition{D}(
        atom::NamedAtom,
        pos::AbstractVector{<:Real},
        occ::Real
    ) where D
        return new(atom, pos, occ)
    end
    function FractionalAtomPosition(
        atom::NamedAtom,
        pos::StaticVector{D,<:Real},
        occ::Real
    ) where D
        return new{D}(atom, pos, occ)
    end
end

"""
    CartesianAtomPosition{D}

Describes an absolute atomic position. The coordinate in the `pos` field is assumed to be given in
angstroms.

Occupancy information is provided in the `occ` field. Note that no checking is done to ensure that
the occupancy is a reasonable value.
"""
struct CartesianAtomPosition{D} <: AbstractAtomPosition{D}
    atom::NamedAtom
    pos::SVector{D,Float64}
    occ::Float64
    function CartesianAtomPosition{D}(
        atom::NamedAtom,
        pos::AbstractVector{<:Real},
        occ::Real
    ) where D
        return new(atom, pos, occ)
    end
    function CartesianAtomPosition(
        atom::NamedAtom,
        pos::StaticVector{D,<:Real},
        occ::Real
    ) where D
        return new{D}(atom, pos, occ)
    end
end

"""
    (T::Type{<:AbstractAtomPosition})(
        [name::AbstractString],
        num::Integer,
        pos::AbstractVector{<:Real},
        occ::Real = 1
    ) -> T

Generates an atomic position from name, atomic number, coordinate, and occupancy data. If a name is
not provided, it's selected automatically from the atomic number.
"""
function (::Type{T})(
    name::AbstractString,
    num::Integer,
    pos::AbstractVector{<:Real},
    occ::Real = 1
) where {T<:AbstractAtomPosition}
    return T(NamedAtom(name, num), pos, occ)
end

function (::Type{T})(
    num::Integer,
    pos::AbstractVector{<:Real},
    occ::Real = 1
) where {T<:AbstractAtomPosition}
    return T(NamedAtom(num), pos, occ)
end

function (::Type{T})(
    name::AbstractString,
    pos::AbstractVector{<:Real},
    occ::Real = 1
) where {T<:AbstractAtomPosition}
    return T(NamedAtom(name), pos, occ)
end

NamedAtom(p::AbstractAtomPosition) = p.atom
displacement(p::AbstractAtomPosition) = p.pos
occupancy(p::AbstractAtomPosition) = p.occ

name(p::AbstractAtomPosition) = name(p.atom)
atomic_number(p::AbstractAtomPosition) = atomic_number(p.atom)
isdummy(p::AbstractAtomPosition) = isdummy(p.atom)

"""
    isapprox(a::AbstractAtomPosition, b::AbstractAtomPosition; atol = sqrt(eps(Float64)))...)

Checks whether two atomic sites are approximately equal to one another. The function returns `true`
if the atomic numbers of the atoms are the same, and the coordinates of the atoms differ by no more
than `atol`. 
"""
function Base.isapprox(a::T, b::T; kwargs...) where T<:AbstractAtomPosition
    a.atom.num == b.atom.num || return false
    return all(isapprox.(a.pos, b.pos; atol=sqrt(eps(Float64)), kwargs...))
end

"""
    CartesianAtomPosition(b::RealBasis{D}, p::FractionalAtomPosition{D})

Converts a fractional atom position to a Cartesian atom position using the supplied basis.
"""
function CartesianAtomPosition(b::RealBasis{D}, p::FractionalAtomPosition{D}) where D
    return CartesianAtomPosition(p.atom, b * p.pos, p.occ)
end

"""
    FractionalAtomPosition(b::RealBasis{D}, p::CartesianAtomPosition{D})

Converts a Cartesian atom position to a fractional atom position using the supplied basis.
"""
function FractionalAtomPosition(b::RealBasis{D}, p::CartesianAtomPosition{D}) where D
    return FractionalAtomPosition(p.atom, b \ p.pos, p.occ)
end

"""
    distance(a1::CartesianAtomPosition, a2::CartesianAtomPosition) -> Float64    
    distance(b::AbstractBasis, a1::FractionalAtomPosition, a2::FractionalAtomPosition) -> Float64

Calculates the distance between two `FractionalAtomPosition` objects in the same basis `b`.
"""
function distance(a1::CartesianAtomPosition, a2::CartesianAtomPosition)
    return norm(displacment(a1) - displacement(a2))
end

function distance(b::AbstractBasis, a1::FractionalAtomPosition, a2::FractionalAtomPosition)
    return norm(RealBasis(b) * (displacement(a1) - displacement(a2)))
end

"""
    deduplicate(l::AbstractVector{T<:AbstractAtomPosition}; atol=sqrt(eps(Float64))) -> Vector{T}
    deduplicate(l::AbstractAtomList; atol=sqrt(eps(Float64))) -> <:AbstractAtomList

Removes atoms that are duplicates or close to being duplicates. In order to be considered
duplicates, the atoms must have both the atomic number, and their coordinates must be approximately
equal (to within a total distance of `sqrt(eps(Float64))`).
"""
function deduplicate(l::AbstractVector{<:AbstractAtomPosition}; atol=sqrt(eps(Float64)))
    # Store the indices of the atoms to keep in a vector
    kept_inds = collect(eachindex(l))
    # Loop through each pair of atoms
    for m in eachindex(l)
        # Only check the remaining atoms if we haven't already excluded this atom
        iszero(m) || for n in m+1:lastindex(l)
            if isapprox(l[m], l[n]; atol)
                # When a duplicate is found, set the value we're checking to zero
                kept_inds[n] = 0
                @debug string(
                    "Removing atom $n:\n",
                    "Atom $m: ", repr("text/plain", l[m]),
                    "Atom $n: ", repr("text/plain", l[n])
                )
            end
        end
    end
    # Filter kept_inds of zeros to get the list of kept atoms
    return l[filter(!iszero, kept_inds)]
end

"""
    move_into_cell(l::AbstractVector{T<:FractionalAtomPosition}; atol=sqrt(eps(Float64)))
        -> Vector{T}
    move_into_cell(l::PeriodicAtomList; atol=sqrt(eps(Float64))) -> PeriodicAtomList

Moves atoms that may exist outside of the bounds of a unit cell (meaning that their fractional
coordinates are not between 0 and 1) into the unit cell.
"""
function move_into_cell(l::AbstractVector{<:FractionalAtomPosition}; atol=sqrt(eps(Float64)))
    return deduplicate(
        [FractionalAtomPosition(NamedAtom(a), mod.(displacement(a), 1), occupancy(a)) for a in l]
    )
end

#---Lists of atoms--------------------------------------------------------------------------------#
"""
    PeriodicAtomList{D}

Contains a list of `FractionalAtomPosition` objects with an associated basis, corresponding to atoms
in a system with periodicity.
"""
struct PeriodicAtomList{D} <: AbstractAtomList{D}
    basis::RealBasis{D}
    atoms::Vector{FractionalAtomPosition{D}}
    function PeriodicAtomList(
        b::AbstractBasis{D},
        l::AbstractVector{FractionalAtomPosition{D}}
    ) where D
        return new{D}(b, deduplicate(l))
    end
end

basis(fl::PeriodicAtomList) = fl.basis

"""
    AtomList{D}

Contains a list of `CartesianAtomPosition` objects, corresponding to atoms in free space without
boundary conditions.
"""
struct AtomList{D} <: AbstractAtomList{D}
    atoms::Vector{CartesianAtomPosition{D}}
    AtomList(l::AbstractVector{CartesianAtomPosition{D}}) where D = new{D}(deduplicate(l))
end

function Base.:(==)(l1::AbstractAtomList, l2::AbstractAtomList) 
    return (l1.basis === l2.basis && l1.atoms == l2.atoms)
end

Base.isempty(l::AbstractAtomList) = isempty(l.atoms)
Base.length(l::AbstractAtomList) = length(l.atoms)
Base.size(l::AbstractAtomList) = size(l.atoms)

Base.getindex(l::AbstractAtomList, ind...) = getindex(l.atoms, ind...)
Base.setindex!(l::AbstractAtomList, x, ind...) = setindex!(l.atoms, x, ind...)
Base.keys(l::AbstractAtomList) = keys(l.atoms)

Base.iterate(l::AbstractAtomList, i::Integer = 1) = iterate(l.atoms, i)

function Base.push!(l::AtomList{D}, a::CartesianAtomPosition{D}) where D
    push!(l.atoms, a)
    return l
end

function Base.push!(l::PeriodicAtomList{D}, a::FractionalAtomPosition{D}) where D
    push!(l.atoms, a)
    return l
end

function Base.push!(l::PeriodicAtomList{D}, a::CartesianAtomPosition{D}) where D
    push!(l.atoms, FractionalAtomPosition(basis(l), a))
    return l
end

Base.pop!(l::AbstractAtomList) = pop!(l.atoms)

Base.filter(f, l::AtomList) = AtomList(filter(f, l.atoms))
Base.filter(f, l::PeriodicAtomList) = PeriodicAtomList(basis(l), filter(f, l.atoms))

function Base.sort!(l::AbstractAtomList; kwargs...)
    sort!(l.atoms; by=NamedAtom, kwargs...)
    return l
end

Base.sort(l::AbstractAtomList; kwargs...) = sort!(deepcopy(l); kwargs...)

"""
    AtomList(l::PeriodicAtomList)

Converts a `PeriodicAtomList` to a list of Cartesian coordinates.
"""
AtomList(l::PeriodicAtomList) = AtomList(map(x -> CartesianAtomPosition(basis(l), x), l))

"""
    PeriodicAtomList(b::RealBasis{D}, l::AbstractVector{CartesianAtomPosition{D}})
    PeriodicAtomList(b::RealBasis{D}, l::AtomList{D})

Uses the supplied basis vectors to convert Cartesian atomic positions in an `AtomList` to fractional
positions with an associated basis.
"""
function PeriodicAtomList(b::AbstractBasis{D}, l::AbstractVector{CartesianAtomPosition{D}}) where D
    return PeriodicAtomList(b, map(x -> FractionalAtomPosition(b,x), l))
end

PeriodicAtomList(b::AbstractBasis{D}, l::AtomList{D}) where D = PeriodicAtomList(b, l.atoms)

deduplicate(l::AtomList; kw...) = AtomList(deduplicate(l.atoms; kw...))
deduplicate(l::PeriodicAtomList; kw...) = PeriodicAtomList(basis(l), deduplicate(l.atoms; kw...))

function move_into_cell(l::PeriodicAtomList; kw...)
    return PeriodicAtomList(basis(l), move_into_cell(l.atoms; kw...))
end

"""
    supercell(l::PeriodicAtomList, M) -> PeriodicAtomList

Creates a new `AtomList` with the basis vectors of a supercell generated by transforming the basis 
vectors of the space by `M`, which may be an integer matrix, an integer vector which is treated as
a diagonal matrix, or a plain integer, which performs a uniform scaling. This function will also
generate new atomic positions to fill the cell.

The function performs this transformation by calculating the Smith normal form of the transformation
matrix. This matrix provides the integer scaling factors needed to stretch the supercell, and the
left unimodular factor is then used to perform the final transformation.
"""
function supercell(l::PeriodicAtomList{D}, M::AbstractMatrix{<:Integer}) where D
    # Convert the provided basis vectors to the supercell basis
    # Conversion to upper triangular form should make this easier to work with
    scb = triangularize(basis(l), M)
    # Get the Smith normal form of the transformation matrix
    (S,U) = snf(M)
    # Generate all the sites in the supercell where new atoms have to be placed
    newpts = vec([SVector(v.I .- 1) for v in CartesianIndices(Tuple(1:s for s in diag(S)))])
    # Shift everything over for each new point
    sclist = [
        FractionalAtomPosition(
            NamedAtom(atom),
            mod.(SMatrix{D,D}(M)\(displacement(atom) + U\d), 1),
            occupancy(atom)
        )
        for atom in move_into_cell(l), d in newpts
    ]
    return PeriodicAtomList(scb, vec(sclist))
end

supercell(l::PeriodicAtomList{D}, x) where D = supercell(l, convert_to_transform(x, Val{D}()))

"""
    atomtypes(l::AbstractAtomList; dummy=false) -> Vector{NamedAtom}

Returns all unique `NamedAtom` types found in an `AbstractAtomList`. This vector is sorted by atomic
number.

The `dummy` keyword controls whether dummy atoms are counted as a separate atom type (`false` by
default).

To obtain a list of all unique atom names or atomic numbers, use `name.(atomtypes(l))` or
`num.(atomtypes(l))`.
"""
function atomtypes(l::AbstractAtomList; dummy::Bool=false)
    return sort!(unique([NamedAtom(a) for a in (dummy ? l : filter(!isdummy, l))]))
end

"""
    atomcounts(l::AbstractAtomList; dummy=false, names=false) -> Vector{Pair{NamedAtom,Int}}

Returns pairs of atoms and the number of atoms in the `AtomList` with that atomic number.

The `dummy` keyword controls whether dummy atoms are counted as a separate atom type (`false` by
default).

The `use_names` keyword determines whether the atoms counted separately based on atom names. By
default, this is equal to `dummy`, so names are only factored in if dummy atoms are counted.
"""
function atomcounts(l::AbstractAtomList; dummy::Bool=false, use_names::Bool = dummy)
    types = atomtypes(l; dummy)
    types = use_names ? types : unique!(reset_name.(types))
    return [Pair(n, count(a -> NamedAtom(a) == n, l)) for n in types]
end

"""
    natomtypes(l::AbstractAtomList; dummy=false) -> Int

Returns the number of types of atoms in an `AbstractAtomList`.

The `dummy` keyword controls whether dummy atoms are counted as a separate atom type (`false` by
default).

The `use_names` keyword determines whether the atoms counted separately based on atom names. By
default, this is equal to `dummy`, so names are only factored in if dummy atoms are counted.
"""
function natomtypes(l::AbstractAtomList; dummy::Bool=false, use_names::Bool = dummy)
    types = atomtypes(l; dummy)
    return length(use_names ? types : unique!(reset_name.(types)))
end
