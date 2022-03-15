"""
    RealSpaceDataGrid{D,T} <: AbstractRealSpaceData{D}

A data grid defined in real space, containing data of type T.
"""
struct RealSpaceDataGrid{D,T} <: AbstractRealSpaceData{D}
    # Basis vectors defining the lattice
    latt::SMatrix{D,D,Float64}
    # Shift of the origin from the lattice
    orig::SVector{D,Float64}
    # The actual data grid
    grid::Array{T,D}
    # Inner constructor
    function RealSpaceDataGrid{D,T}(
        latt::AbstractMatrix{<:Real},
        orig::AbstractVector{<:Real},
        grid::AbstractArray{T,D}
    ) where {D,T}
        # Still gotta pass the lattice sanity checks
        lattice_sanity_check(latt)
        return new(latt, orig, grid)
    end
end

"""
    RealSpaceDataGrid{D,T}(latt::AbstractLattice{D}, grid::AbstractArray{T,D}; prim=true)

Creates a real space data grid using lattice information from an `AbstractLattice`. By default,
data is assumed to be given in terms of the primitive lattice (as is usually the case for
computational data).
"""
function RealSpaceDataGrid{D,T}(
    latt::AbstractLattice{D},
    orig::AbstractVector{<:Real},
    grid::AbstractArray{T,D};
    prim=true
) where {D,T}
    # Conversion for safety
    l = RealLattice{D}(latt)
    if prim
        return RealSpaceDataGrid{D,T}(prim(l), orig, grid)
    else
        return RealSpaceDataGrid{D,T}(conv(l), orig, grid)
    end
end

# Data in RealSpaceDataGrids can now be indexed
Base.getindex(g::RealSpaceDataGrid{D,T} where {D,T}, inds...) = getindex(g.grid, inds...)

"""
    basis(g::RealSpaceDataGrid{D,T}) -> SMatrix{D,D,Float64}

Gets the basis vectors of a `RealSpaceDataGrid`.
"""
basis(g::RealSpaceDataGrid) = g.latt
grid(g::RealSpaceDataGrid) = g.grid
# Size of the data grid in entries per dimension
# TODO: should we overload Base.size() as well?
gridsize(g::RealSpaceDataGrid{D,T} where {D,T}) = size(g.grid)
#Base.size(g::RealSpaceDataGrid{D,T} where {D,T}) = gridsize(g)

# TODO: Make mathematical operations on RealSpaceDataGrids easy

"""
    KPointGrid{D} <: AbstractKPoints{D}

Contains a grid used to generate k-points during a calculation.

The grid itself is given as an `SMatrix{D,D,Int}`, and can be interpreted as a set of `D` vectors 
given in terms of the primitive basis. These vectors can alternatively be used to construct a 
supercell.

The shift of the k-point mesh off Γ is given as an `SVector{D,Float64}`.
"""
struct KPointGrid{D} <: AbstractKPoints{D}
    grid::SMatrix{D,D,Int}
    orig::SVector{D,Float64}
    function KPointGrid{D}(grid::AbstractMatrix{<:Integer}, orig::AbstractVector{<:Real}) where D
        # only allow positive values in the grid matrix
        @assert all(x -> x > 0, grid) "negative values are disallowed in the grid matrix"
        # Keep shift inside the Brillouin zone
        orig = orig -  round.(orig)
        return new(grid, orig)
    end
end

"""
    KPointList{D} <: AbstractKPoints{D}

Contains a list of k-points. This is useful for describing an ordered list of k-points that are 
not associated with a mesh - for instance, those used in band structure calculations. This can
also be used to store lists of k-points that are generated by a `KPointGrid{D}`.
"""
struct KPointList{D} <: AbstractKPoints{D}
    list::Vector{SVector{D,Float64}}
end

# Index like the vector it is internally
Base.getindex(k::KPointList{D}, i) where D = k.list[i]

"""
    nkpt(k::KPointList{D}) -> Int

Gets the number of k-points in a `KPointList`.
"""
nkpt(k::KPointList{D}) where D = length(k.list)

#= TODO: figure out how to get a k-point list from a grid
function KPointList{D}(k::KPointGrid)

end
=#

"""
    BandAtKPoint

Stores information about a band's energy and its occupancy at a specific k-point.
"""
struct BandAtKPoint
    # Energies
    e::Vector{Float64}
    # Occupancy
    occ::Vector{Float64}
    function BandAtKPoint(e::AbstractVector{<:Real}, occ::AbstractVector{<:Real})
        @assert length(e) == length(occ) "Size of energy and occupancy arrays do not match."
        return new(e, occ)
    end
end

"""
    BandAtKPoint(eocc::AbstractVector{NTuple{2,<:Real}})

Constructs a new `BandAtKPoint` from a vector containing tuples of energy and occupancy data
(in that order).
"""
function BandAtKPoint(eocc::AbstractVector{<:NTuple{2,<:Real}})
    return BandAtKPoint([x[1] for x in eocc], [x[2] for x in eocc])
end

# Access any pair of energy and occupancy with indexing
function Base.getindex(b::BandAtKPoint, inds...)
    return (b.e[inds...], b.occ[inds...])
end

"""
    nband(b::BandAtKPoint) -> Int

Returns the number of bands associated with a k-point.
"""
nband(b::BandAtKPoint) = length(b.e)

"""
    BandStructure{D}

Stores information about an electronic band structure, including the list of k-points used to
generate the data (as a `KPointList{D}`)and the band information at every k-point (as a 
`Vector{BandAtKPoint}`).
"""
struct BandStructure{D} <: AbstractReciprocalSpaceData{D}
    # k-points for which band data is defined
    kpts::KPointList{D}
    # Set of energy and occupancy data
    bands::Vector{BandAtKPoint}
    function BandStructure{D}(kpts::KPointList{D}, bands::Vector{BandAtKPoint}) where D
        @assert nkpt(kpts) == length(bands) "Incorrect number of k-points or band datasets."
        @assert _allsame(length(bands)) "Number of bands is inconsistent."
        return new(kpts, bands)
    end
end

"""
    BandStructure{D}(kpts::AbstractKPoints{D}, bands::AbstractVector{<:BandAtKPoint}) where D

Generates a new band structure from k-point information and a vector containing band information
at each k-point.
"""
function BandStructure{D}(kpts::AbstractKPoints{D}, bands::AbstractVector{<:BandAtKPoint}) where D
    return BandStructure{D}(kpts, bands)
end

# Get the pair of a k-point and associated band data
function Base.getindex(b::BandStructure{D}, inds...) where D
    return (b.kpts[inds...], b.bands[inds...])
end

nkpt(b::BandStructure{D}) where D = nkpt(b.kpts)
nband(b::BandStructure{D}) where D = nband(b.bands[1])

"""
    HKLData{D,T} <: AbstractReciprocalSpaceData{D}

Stores information associated with specific sets of reciprocal lattice vectors. Data can be
accessed and modified using regular indexing, where indices may be negative.
"""
struct HKLData{D,T} <: AbstractReciprocalSpaceData{D}
    # the actual data
    data::Array{T,D}
    # the bounds in each dimension
    # mutable since the dimensions of Array{D,T} can be changed, in principle
    bounds::MVector{D,UnitRange{Int}}
    function HKLData{D,T}(
        data::AbstractArray{T,D},
        bounds::AbstractVector{<:AbstractRange{<:Integer}}
    ) where {D,T}
        # The size of the array should match the bounds given
        # For instance, a HKLData with bounds  [-10:10, -10:10, -10:10] should Be
        # [21, 21, 21]
        @assert [s for s in size(data)] == [length(r) for r in bounds] "Array size \
        incompatible with bounds."
        return new(data, bounds)
    end
end

# Needed because HKLData will nearly always have unexpected indices
Base.has_offset_axes(g::HKLData) = true

"""
    HKLData_boundscheck(g::HKLData{D,T}, inds::Vararg{<:Integer,D})

Checks that array indices used to access data in an `HKLData` are valid.
"""
function shiftbounds(g::HKLData{D,T}, inds) where {D,T}
    # Check that all the indices are in bounds
    if !mapreduce((i,r) -> i in r, &, inds, g.bounds)
        # TODO: does this produce a reasonable error message with correct bounds?
        throw(BoundsError(g, inds))
    end
    # Adjust the indices to match the array
    # Subtract the minimum index then add 1
    # So if the range is -10:10, an index of 0 should be 0 - -10 + 1 = 11
    i = inds .- minimum.(g.bounds) .+ 1
    return i
end

# HKLData now supports indexing by Miller index
function Base.getindex(g::HKLData{D,T}, inds...) where {D,T}
    i = shiftbounds(g, inds)
    return g.data[i...]
end

function Base.setindex!(g::HKLData{D,T}, x::T, inds...) where {D,T}
    i = shiftbounds(g, inds)
    g.data[i...] = x
end

function HKLData(
    a::AbstractArray{T},
    bounds::Vararg{AbstractUnitRange{<:Integer},D}
) where {D,T}
    return HKLData{D,T}(a, bounds)
end

function Base.zeros(
    ::Type{HKLData{D,T}},
    bounds::Vararg{AbstractUnitRange{<:Integer}, D}
) where {D,T}
    data = zeros(T, length.(bounds))
    return HKLData{D,T}(data, MVector{D,UnitRange{Int}}(bounds))
end

"""
    ReciprocalWavefunction{D,T<:Real} <: AbstractReciprocalSpaceData{D}

Contains a wavefunction stored by k-points and bands in a planewave basis. Used to store data in
VASP WAVECAR files.

For every k-point, there is an associated number of bands (the same number for every k-point). Band
information is stored as the energy of the band and its occupancy.

Every band has associated data containing coefficients of the constituent planewaves stored in a 
`HKLData{D,Complex{T}}`. Unlike most data structures provided by this package, the type of
complex number used does not default to `Float64`: wavefunction data is often supplied as a 
`Complex{Float32}` to reduce the size of the data.
"""
struct ReciprocalWavefunction{D,T<:Real} <: AbstractReciprocalSpaceData{D}
    # Reciprocal lattice on which the k-points are defined
    rlatt::SMatrix{D,D,Float64}
    # Band information (energy and occupation) at each k-point
    bands::BandStructure{D}
    # Planewave coefficients
    # Vector (size nkpt) of Vectors (size nband) of HKLData
    waves::Vector{Vector{HKLData{D,Complex{T}}}}
    function ReciprocalWavefunction{D,T}(
        rlatt::AbstractMatrix{<:Real},
        bands::BandStructure{D},
        waves::AbstractVector{<:AbstractVector{HKLData{D,Complex{T}}}}
    ) where {D,T<:Real}
        # Number of k-points should equal the length of the waves vector
        @assert length(waves) == nkpt(bands) "Number of k-points and number of planewave data \
        entries do not match."
        # Number of band entries per k-point should be the same
        lw = length.(waves)
        @assert _allsame(lw) "The number of bands per k-point is inconsistent."
        # There should be the same number of planewave sets as there are bands per k-point``
        @assert first(lw) == nband(bands) "Number of bands and number of planewave data \
        entries do not match."
        return new(rlatt, bands, waves)
    end
end

function ReciprocalWavefunction{D,T}(
    latt::AbstractLattice{D},
    bands::BandStructure{D},
    waves::AbstractVector{<:AbstractVector{HKLData{D,Complex{T}}}};
    prim = true
) where {D,T<:Real}
    M = pick(latt, prim)
    return ReciprocalWavefunction{D,T}(M, bands, waves)
end

# Getting indices should pull from the waves struct
# This pulls first from k-points, then from bands, then from HKL data
getindex(wf::ReciprocalWavefunction{D,T}, inds...) where {D,T} = wf.waves[inds...]

nkpt(wf::ReciprocalWavefunction{D,T}) where {D,T} = nkpt(wf.bands)
nband(wf::ReciprocalWavefunction{D,T}) where {D,T} = nband(wf.bands)

"""
    DensityOfStates

Contains the total density of states information.
"""
struct DensityOfStates <: AbstractDensityOfStates
    # Fermi energy
    fermi::Float64
    # Energy at each point
    energy::Vector{Float64}
    # Density of states
    dos::Vector{Float64}
    # Integrated density of states
    int::Vector{Float64}
    function DensityOfStates(
        fermi::Real,
        energy::AbstractVector{<:Real},
        dos::AbstractVector{<:Real},
        int::AbstractVector{<:Real}
    )
        @assert size(energy) == size(dos) == size(int) string(
            "Sizes of input vectors do not match!"
        )
        return new(fermi, energy, dos, int)
    end
end

"""
    fermi(d::AbstractDensityOfStates) -> Float64

Gets the Fermi energy from DOS data. There are no guarantees on the unit of energy used!
"""
fermi(d::AbstractDensityOfStates) = d.fermi

"""
    ProjectedDensityOfStates

Contains projected density of states information.
"""
struct ProjectedDensityOfStates
    fermi::Float64
    energy::Vector{Float64}
    dos::Matrix{Float64}
end