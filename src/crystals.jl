"""
    Crystal{D} <: AbstractCrystal{D}

A crystal structure in `D` dimensions. Contains information about the lattice, space group, and
atoms.

Atoms can be specified in terms of the generating set, which can be used to fill a cell given 
space group information, or in terms of explicitly given atomic positions.

If `sgno` is set to 0, it will be assumed that coordinates are given in Cartesian coordinates in 
units of angstroms. Otherwise, it will be assumed that coordinates are given in terms of the unit
cell basis.
"""
struct Crystal{D} <: AbstractCrystal{D}
    # Lattice information (real space)
    latt::RealLattice{D}
    # Space group number - set to 0 if non-periodic/unknown
    sgno::Int
    # Space group origin (location of inversion center, if present - defaults to zeros)
    orig::SVector{D,Float64}
    # Positions of generating atoms (needed to fill the whole structure given the space group)
    gen::AtomList{D}
    # Positions of atoms explicitly generated (used to generate templates, XYZs, etc.)
    pos::AtomList{D}
    # Inner constructor
    function Crystal(
        latt::AbstractLattice{D},
        sgno::Integer,
        orig::AbstractVector{<:Real},
        gen::AtomList{D},
        pos::AtomList{D}
    ) where D
        # TODO: include some validation
        # If the atom list doesn't have a basis defined (assuming Cartesian coordinates)
        # generate a new AtomList 
        if basis(gen) == zeros(BasisVectors{D})
            # Do we want to use the conventional basis vectors all the time?
            # I think it's a good temporary choice, just because if a dataset is being loaded in
            # from a computation, conv(latt) should generally match the primitive cell
            gen = AtomList{D}(conv(latt), reduce_coords(conv(latt), gen, incell=true))
        end
        return new{D}(RealLattice(latt), sgno, orig, gen, pos)
    end
end

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
Base.propertynames(::CrystalWithDatasets) = tuple(
    fieldnames(CrystalWithDatasets)...,
    fieldnames(Crystal)...,
    fieldnames(Dict)...
)

function Base.getproperty(xtal::CrystalWithDatasets, f::Symbol)
    # Core fields
    if f in fieldnames(CrystalWithDatasets)
        return getfield(xtal, f)
    # Fields in the Crystal{D}
    elseif f in fieldnames(Crystal)
        return getproperty(getfield(xtal, :xtal), f)
    # Fields in the Dict{K,V}
    elseif f in fieldnames(Dict)
        return getproperty(getfield(xtal, :data), f)
    else
        error("type CrystalWithDatasets has no field ", string(f))
    end
end

# Allow for getting datasets by key; no need to reach into the Dict
# TODO: figure out how autocompletion works and how to enable it here
# UPDATE: seems like it's hard-coded into the REPL for `AbstractDict`
# So adding in that functionality doesn't seem to be feasible
function Base.getindex(xtaldata::CrystalWithDatasets{D,K,V}, key::K) where {D,K,V}
    return xtaldata.data[key]
end

# Easy way of pulling just the crystal from a `CrystalWithDatasets{D}`
function Base.convert(::Type{Crystal{D}}, xtaldata::CrystalWithDatasets{D,K,V}) where {D,K,V}
    return xtaldata.xtal
end

Crystal(xtaldata::CrystalWithDatasets{D,K,V}) where {D,K,V} = xtaldata.xtal
data(xtaldata::CrystalWithDatasets{D,K,V}) where {D,K,V} = xtaldata.data

RealLattice(xtal::AbstractCrystal) = xtal.latt
ReciprocalLattice(xtal::AbstractCrystal) = ReciprocalLattice(xtal.latt)

prim(xtal::AbstractCrystal) = prim(RealLattice(xtal))
conv(xtal::AbstractCrystal) = conv(RealLattice(xtal))

basis(xtal::AbstractCrystal; primitive::Bool=false) = primitive ? prim(xtal.latt) : conv(xtal.latt)

volume(xtal::AbstractCrystal; primitive::Bool=false) = volume(xtal.latt; primitive)

atoms(xtal::AbstractCrystal) = xtal.pos

# TODO: fix this so that it gets the right number of atoms regardless of space group
"""
    natom_cell(xtal::Crystal)

Returns the number of atoms in a crystal's unit cell.

This function is currently not aware of space groups or settings, so if the generating set does not
contain all the atoms in the cell, it will return the wrong value.
"""
natom_cell(xtal::AbstractCrystal) = natom(xtal.gen)

"""
    natom_template(xtal::Crystal)

Returns the number of atoms in the crystal template.
"""
natom_template(xtal::AbstractCrystal) = natom(xtal.pos)

atomtypes(xtal::AbstractCrystal) = atomtypes(xtal.gen)
natomtypes(xtal::AbstractCrystal) = natomtypes(xtal.gen)

#=
"""
    all_atoms_in_cell(spgrp::Integer, l::AtomList{D}; onbounds::Bool = false) -> AtomList{D}

Generates all atoms within a unit cell given the space group information and a minimal set of
atomic positions.
"""
function all_atoms_in_cell(spgrp::Integer, l::AtomList{D}; onbounds::Bool = false) where D
    # Space group 0 effectively means no translation symmetry
    # Space group 1 should be the trivial group (just translations): do nothing
    spgrp in (0,1) && return l
    # TODO: implement this function.
    return l
end
=#
