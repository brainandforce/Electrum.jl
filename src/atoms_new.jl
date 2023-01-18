"""
    NamedAtom

Stores information about an atom, which includes a name which may be up to 8 codepoints long, and
the atomic number.

Internally, the name is stored as a `NTuple{8,Char}` in the `name` field to guarantee that the type
is pure bits. However, the `name` property returns a `String`.
"""
struct NamedAtom
    name::NTuple{8,Char}
    num::Int
    function NamedAtom(name::AbstractString, num::Integer)
        codepoints = ntuple(Val{8}()) do i
            i <= length(name) ? name[i] : '\0'
        end
        return new(codepoints, num)
    end
end

function Base.getproperty(atom::NamedAtom, p::Symbol)
    if p == :name
        l = findfirst(isequal('\0'), getfield(atom, name))
        return string(getfield(atom, name)[1:l-1])
    end
    return getfield(atom, p)
end

NamedAtom(num::Integer) = num in 1:118 ? NamedAtom("dummy", num) : NamedAtom(ELEMENTS[num], num)

Base.show(io::IO, atom::NamedAtom) = println(io, "NamedAtom(", atom.name, ", ", atom.num, ")")
Base.isless(a1::NamedAtom, a2::NamedAtom) = isless(a1.num, a2.num) || isless(a1.name, a2.name)

name(a::NamedAtom) = a.name
atomic_number(a::NamedAtom) = a.num
isdummy(a::NamedAtom) = (a.num == 0)

reset_name(a::NamedAtom) = NamedAtom(atomic_number(a))

#---Atomic position data--------------------------------------------------------------------------#
"""
    AbstractAtomPosition{D}

Supertype that describes atomic positions in `D` dimensions, which include name, coordinate, and
occupancy information.
"""
abstract type AbstractAtomPosition{D}
end

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

NamedAtom(p::AbstractAtomPosition) = p.atom
position(p::AbstractAtomPosition) = p.pos
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
    FractionalAtomPosition(b::RealBasis{D}, p::FractionalAtomPosition{D})

Converts a Cartesian atom position to a fractional atom position using the supplied basis.
"""
function FractionalAtomPosition(b::RealBasis{D}, p::FractionalAtomPosition{D}) where D
    return FractionalAtomPosition(p.atom, b \ p.pos, p.occ)
end

#---Lists of atoms--------------------------------------------------------------------------------#
"""
    AbstractAtomList{D}

Supertype for lists of atomic positions in `D` dimensions.
"""
abstract type AbstractAtomList{D}
end

"""
    PeriodicAtomList{D}

Contains a list of `FractionalAtomPosition` objects with an associated basis, corresponding to
atoms in a system with periodicity.
"""
struct PeriodicAtomList{D} <: AbstractAtomList{D}
    basis::RealBasis{D}
    atoms::Vector{FractionalAtomPosition{D}}
    function PeriodicAtomList(
        b::AbstractBasis{D},
        l::AbstractVector{FractionalAtomPosition{D}}
    ) where D
        return new{D}(b,l)
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
    AtomList(l::AbstractVector{CartesianAtomPosition{D}}) where D = new(b,l)
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

Base.filter(f, l::AtomList) = AtomList(filter(f, l.atoms))
Base.filter(f, l::PeriodicAtomList) = PeriodicAtomList(basis(l), filter(f, l.atoms))

Base.sort(l::AtomList; kwargs...) = AtomList(sort(l.atoms); by=NamedAtom, kwargs...)

function Base.sort(l::PeriodicAtomList; kwargs...)
    return PeriodicAtomList(basis(l), sort(l.atoms); by=NamedAtom, kwargs...)
end

function AtomList(l::PeriodicAtomList{D}) where D
    return AtomList(CartesianAtomPosition.(basis(l), l.atoms))
end

"""
    deduplicate(l::AbstractVector{T<:AbstractAtomPosition}; atol=sqrt(eps(Float64))) -> Vector{T}

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

deduplicate(l::AtomList; kw...) = AtomList(deduplicate(l.atoms; kw...))
deduplicate(l::PeriodicAtomList; kw...) = PeriodicAtomList(basis(l), deduplicate(l.atoms; kw...))

"""
    supercell(l::AtomList, M::AbstractMatrix{<:Integer}) -> AtomList

Creates a new `AtomList` with the basis vectors of a supercell generated by transforming the basis 
vectors of the space by `M`. This will also generate new atomic positions to fill the cell.
"""
function supercell(l::PeriodicAtomList{D}, M::AbstractMatrix{<:Integer}) where D
    # Convert the provided basis vectors to the supercell basis
    # Conversion to upper triangular form should make this easier to work with
    scb = triangularize(basis(l), M)
    # Get the Smith normal form of the transformation matrix
    (S,U) = snf(M)
    # Generate all the sites in the supercell where new atoms have to be placed
    newpts = vec([SVector(v.I .- 1) for v in CartesianIndices(Tuple(diag(S)))])
    # Move all positions into the cell; remove duplicates
    dedup = deduplicate(
        [FractionalAtomPosition(name(a), atomicno(a), mod.(position(a), 1)) for a in l]
    )
    # Shift everything over for each new point
    sclist = [
        AtomPosition(
            name(atom),
            atomicno(atom),
            # Keep the new atoms inside the supercell
            # Multiply the coordinate by U
            mod.(SMatrix{D,D}(M)\(position(atom) + U*d), 1)
        )
        for atom in dedup, d in newpts
    ]
    return AtomList(scb, vec(sclist))
end

supercell(l::AbstractAtomList, v::AbstractVector{<:Integer}) = supercell(l, diagm(v))

"""
    atomtypes(l::AbstractAtomList; dummy=false) -> Vector{NamedAtom}

Returns all unique `NamedAtom` types found in an `AbstractAtomList`. This vector is sorted by
atomic number.

The `dummy` keyword controls whether dummy atoms are counted as a separate atom type (`false` by
default).

To obtain a list of all unique atom names or atomic numbers, use `name.(atomtypes(l))` or
`num.(atomtypes(l))`.
"""
function atomtypes(l::AbstractAtomList; dummy::Bool=false)
    return sort!(unique([NamedAtom(a) for a in (dummy ? filter(isdummy, l) : l)]))
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
    types = use_names ? types : unique!(reset_names.(types))
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
    return length(use_names ? types : unique!(reset_names.(types)))
end

export NamedAtom, AbstractAtomPosition, FractionalAtomPosition, CartesianAtomPosition,
       AbstractAtomList, AtomList, PeriodicAtomList
export name, atomic_number, isdummy, position, occupancy, deduplicate, supercell, atomtypes,
       atomcounts, natomtypes
