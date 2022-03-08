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
    function Crystal{D}(
        latt::AbstractLattice{D},
        sgno::Integer,
        orig::AbstractVector{<:Real},
        gen::AtomList{D},
        pos::AtomList{D}
    ) where D
        return new(RealLattice{D}(latt), sgno, orig, gen, pos)
    end
end

# The name `CrystalData{D}` was not used to avoid implying that this is a subtype of 
# `AbstractCrystalData{D}` - which is used only for datasets
"""
    CrystalWithDatasets{D,K,V} <: AbstractCrystal{D}

A pairing of a `Crystal{D}` and a `Dict{K,V}` which allows for access to associated datasets.
"""
struct CrystalWithDatasets{D,K,V} <: AbstractCrystalData{D}
    xtal::Crystal{D}
    data::Dict{K,V}
end

# Allow for getting datasets by key; no need to reach into the Dict
function Base.getindex(xtaldata::CrystalWithDatasets{D,K,V}, key::K) where {D,K,V}
    return xtaldata.data[key]
end

# Easy way of pulling just the crystal from a `CrystalWithDatasets{D}`
Crystal{D}(xtaldata::CrystalWithDatasets{D,K,V}) where {D,K,V} = xtaldata.xtal