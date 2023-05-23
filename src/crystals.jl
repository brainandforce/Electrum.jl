"""
    AbstractCrystal{D}

A crystal structure in `D` dimensions, containing information about the lattice, space group, and
atoms contained within the crystal.
"""
abstract type AbstractCrystal{D}
end

"""
    Crystal{D} <: AbstractCrystal{D}

A crystal structure in `D` dimensions. Contains information about the lattice, space group, and
atoms. This is a mutable data structure.

At minimum, a list of atomic positions (as an `AtomList`) is needed to generate a `Crystal`.
Optionally, space group number and the origin of the space group may be provided.

A transform may also be specified that converts the basis vectors of the `AtomList` to a favored
representation, most often the conventional lattice. If it is not specified, it is filled with an
identity matrix by default. The matrix is right-multiplied with the basis vectors to produce the
favored representation. Because the rows of the transformation naturally correspond to the
operations performed on each constituent basis vector, it may be easier to enter the transform as a
transpose (or equivalently, an adjoint) when entered manually.
"""
mutable struct Crystal{D} <: AbstractCrystal{D}
    # Positions of generating atoms (needed to fill the whole structure given the space group)
    atoms::PeriodicAtomList{D}
    # Space group number - set to 0 if non-periodic/unknown
    sgno::Int
    # Space group origin (location of inversion center, if present - defaults to zeros)
    sgorig::SVector{D,Float64}
    # Transform to generate the conventional or primitive lattice
    # It might be easier to write this as a transpose
    transform::SMatrix{D,D,Int}
    function Crystal(
        # Always required, even if generating an empty list of atoms
        # (because the included basis is necessary)
        atoms::PeriodicAtomList{D},
        # Assume space group is zero if not provided
        sgno::Integer = 0,
        # Assume the origin is [0, 0, 0] if not provided
        sgorig::AbstractVector{<:Real} = zeros(SVector{D,Float64}),
        # Assume there is no conversion unless otherwise specified
        transform = LinearAlgebra.I
    ) where D
        return new{D}(atoms, sgno, sgorig, convert_to_transform(transform, Val{D}()))
    end
end

function Base.:(==)(x1::T, x2::T) where T<:Crystal
    return all(getfield(x1, s) == getfield(x2, s) for s in fieldnames(T))
end

function Base.getproperty(xtal::Crystal, s::Symbol)
    return s === :basis ? basis(getfield(xtal, :atoms)) : getfield(xtal, s)
end

Base.propertynames(::Crystal) = (:basis, fieldnames(Crystal)...)

# The name `CrystalData{D}` was not used to avoid implying that this is a subtype of 
# `AbstractCrystalData{D}` - which is used only for datasets
"""
    CrystalWithDatasets{D,K,V} <: AbstractCrystal{D}

A pairing of a `Crystal{D}` and a `Dict{K,V}` which allows for access to associated datasets.
"""
struct CrystalWithDatasets{D,K,V} <: AbstractCrystal{D}
    xtal::Crystal{D}
    data::Dict{K,V}
end

# Expose the constituent crystal's fields
Base.propertynames(x::CrystalWithDatasets) = (fieldnames(typeof(x))..., propertynames(x.xtal)...)

function Base.getproperty(xtal::CrystalWithDatasets, s::Symbol)
    # Fields in the Crystal{D}
    s in propertynames(getfield(xtal, :xtal)) && return getproperty(getfield(xtal, :xtal), s)
    return getfield(xtal, s)
end

function Base.setproperty!(xtal::CrystalWithDatasets, s::Symbol, x)
    if s in propertynames(getfield(xtal, :xtal))
        setproperty!(getfield(xtal, :xtal), s, x)
    else
        setfield!(xtal, s, x)
    end
end

# Allow for getting datasets by key; no need to reach into the Dict
# TODO: figure out how autocompletion works and how to enable it here
# UPDATE: seems like it's hard-coded into the REPL for `AbstractDict`
# So adding in that functionality doesn't seem to be feasible
Base.getindex(xtaldata::CrystalWithDatasets, key) = xtaldata.data[key]

# Easy way of pulling just the crystal from a `CrystalWithDatasets{D}`
Base.convert(T::Type{<:Crystal}, xtaldata::CrystalWithDatasets) = xtaldata.xtal::T

Crystal(xtaldata::CrystalWithDatasets) = xtaldata.xtal
data(xtaldata::CrystalWithDatasets) = xtaldata.data

"""
    generators(xtal::AbstractCrystal{D}) -> PeriodicAtomList{D}

Returns the list of generating atomic positions associated with a `Crystal` or
`CrystalWithDatasets`.

Note that this does not convert the input to a `PeriodicAtomList` with all atomic positions; only
the minimal set that's needed to generate the atoms given the space group symmetry. To enumerate
all of the atoms, use `convert(PeriodicAtomList, xtal)` or `PeriodicAtomList(xtal)`.
"""
generators(xtal::AbstractCrystal) = xtal.atoms

function Base.convert(T::Type{<:PeriodicAtomList}, xtal::AbstractCrystal)
    # TODO: when space groups are implemented, explicitly generate those positions
    return supercell(xtal.atoms, xtal.transform)::T
end

"""
    PeriodicAtomList(xtal::AbstractCrystal{D}) -> PeriodicAtomList{D}

Converts a `Crystal` or `CrystalWithDatasets` to a list of all generated atomic positions given
the space group and the associated transformation matrix.

To recover only the generating set of atoms, use `generators(xtal)`.
"""
PeriodicAtomList(xtal::AbstractCrystal) = convert(PeriodicAtomList, xtal)
AtomList(xtal::AbstractCrystal) = AtomList(PeriodicAtomList(xtal))

# TODO: Generate the full atom list when symmetry operations are implemented.
Base.length(xtal::AbstractCrystal) = length(xtal.atoms)

atomtypes(xtal::AbstractCrystal; kwargs...) = atomtypes(PeriodicAtomList(xtal); kwargs...)
atomcounts(xtal::AbstractCrystal; kwargs...) = atomcounts(PeriodicAtomList(xtal); kwargs...)
natomtypes(xtal::AbstractCrystal; kwargs...) = natomtypes(PeriodicAtomList(xtal); kwargs...)

"""
    get_transform(xtal::AbstractCrystal{D}) -> SMatrix{D,D,Int}

Returns the transformation matrix associated with an `AbstractCrystal`.
"""
get_transform(xtal::AbstractCrystal) = xtal.transform

"""
    set_transform!(xtal::AbstractCrystal, M) -> AbstractCrystal

Sets the transform supplied with an `AbstractCrystal`. The transform can be an integer matrix,
vector, scalar, or `UniformScaling`, which is converted to an `SMatrix{D,D,Int}` when stored.

The function returns the modified input for convenience.
"""
function set_transform!(xtal::AbstractCrystal, M)
    setproperty!(xtal, :transform, convert_to_transform(M))
    return xtal
end
