#---General functions------------------------------------------------------------------------------#
"""
    nspin(x) -> Int

Returns the number of spin states associated with a dataset. For 3D data, this is usually 1 or 2,
depending on whether a restricted or unrestricted calculation was run.
"""
function nspin end

"""
    nband(x) -> int

Returns the number of bands associated with a dataset.
"""
function nband end

"""
    ecut(x) -> Union{Real,Missing}

Returns the energy cutoff (in Hartree) associated with a dataset. By default, this returns the value
of `x.ecut`, unless this value is `NaN`, in which case it returns `missing`.
"""
ecut(x) = ifelse(isnan(x.ecut), missing, x.ecut)

"""
    fermi(x) -> Union{Real,Missing}

Returns the Fermi energy (in Hartree) associated with a dataset. By default, this returns the value
of `x.fermi`, unless this value is `NaN`, in which case it returns `missing`.
"""
fermi(x) = ifelse(isnan(x.fermi), missing, x.fermi)

#---Band structures--------------------------------------------------------------------------------#
"""
    Electrum._assert_fermi_energy(data::AbstractArray{<:EnergyOccupancy}, fermi::Real) -> Nothing

Asserts that the Fermi energy lies between the minimum and maximum energy of an `EnergyOccupancy`
array.
"""
function _assert_fermi_energy(data::AbstractArray{<:EnergyOccupancy}, fermi::Real)
    isnan(fermi) && return nothing
    @assert min_energy(data) <= fermi "Fermi energy is below the minimum energy in the dataset."
    @assert fermi <= max_energy(data) "Fermi energy is above the maximum energy in the dataset."
end

"""
    BandStructure{D,T<:Real} <: AbstractArray{EnergyOccupancy{T},3}

Stores band structure information, consisting of `EnergyOccupancy` data in a three-dimensional
row-major matrix of spin, k-point, and band indices (in this order), along with the associated
k-points, energy cutoff, and Fermi energy.

# Indexing

By convention, electronic data is usually indexed by spin, then k-point, then band, but it is more
common to deal with wavefunctions of specific spins separately. Therefore, row-major indexing is
used so that data for spin-up and spin-down states are stored continguously.

# Energy cutoff and Fermi energy

Defining fields named `ecut` and `fermi` containing the energy cutoff and Fermi energy respectively
will automatically define `ecut(::ElectronicEnergyData)` and `fermi(::ElectronicEnergyData)` for 
those types.

## Missing data conventions

`NaN` entries for the `ecut` and `fermi` fields are considered equivalent to `missing` data, and the
accessor functions determine this automatically.
"""
struct BandStructure{D,T} <: AbstractArray{EnergyOccupancy{T},3}
    data::EnergyOccupancy{T,3}
    kpoints::KPointMesh{D}
    ecut::T
    fermi::T
    function BandStructure(
        data::AbstractArray{EnergyOccupancy{S}},
        kpoints::AbstractVector{<:KPoint{D}},
        ecut::Real,
        fermi::Real
    ) where {D,S}
        T = promote_type(S, typeof(ecut), typeof(fermi))
        _assert_fermi_energy(data, fermi)
        @assert ecut > min_energy(data) "Energy cutoff ($ecut Ha) is unreasonably low."
        @assert size(data, 2) == length(kpoints) string(
            "Expected $(size(data, 2)) k-points, got $(length(kpoints))"
        )
        return new{D,T}(data, kpoints, ecut, fermi)
    end
end

Base.size(b::ElectronicEnergyData) = reverse(size(b.data))
Base.getindex(b::ElectronicEnergyData, i...) = getindex(b.data, reverse(i)...)

nspin(b::BandStructure) = size(b, 1)
nband(b::BandStructure) = size(b, 3)

"""
    nkpt(b::BandStructure) -> Int

Returns the number of k-points associated with `b`.
"""
function nkpt(b::BandStructure)
    @assert size(b, 2) == length(b.kpoints) "Mismatch in expected number of k-points"
    return size(b, 2)
end

KPointMesh(b::BandStructure) = b.kpoints

#---Density of states------------------------------------------------------------------------------#
"""
    DensityOfStates{T<:Real} <: AbstractMatrix{StateDensity{T}}

Stores density of states information in a matrix indexed by spin and energy level, as well as the
associated energy cutoff and Fermi energy, if provided.

# Indexing

As with `BandStructure{T}`, the storage mode of `DensityOfStates{T}` is row-major, meaning that the
density of states data associated with each spin is stored contiguously, but the spin index comes
first, then energy index.
"""
struct DensityOfStates{T} <: AbstractMatrix{StateDensity{T}}
    states::Matrix{StateDensity{T}}
    ecut::T
    fermi::T
end

Base.size(d::DensityOfStates) = reverse(size(d.states))
Base.getindex(d::DensityOfStates, i...) = getindex(d.states, reverse(i)...)

"""
    ProjectedDensityOfStates{T<:Real} <: Vector{DensityOfStates{T}}

Stores sets of `DensityOfStates{T}` information that has been projected by atom, orbital, or other
method.

Along with an array of all `StateDensity{T}` information, a set of atom labels stored as 
`InlineStrings.InlineString15` objects (identical to how `NamedAtom` stores its atom labels).

# Indexing

A `ProjectedDensityOfStates{T}` object stores each set of constituent `DensityOfStates{T}` data
contiguously, and indexing is row based, meaning that the first index argument is the atom number or
label, the second is the spin index, and the remaining indices correspond to the energy levels.

The total density of states, if present, resides at integer index 0, or string index `total`.

# Total density of states

By convention, the total density of states is stored at index 0. Whether this data is present is
determined by the the result of `has_total_dos(p::ProjectedDensityOfStates)`.
"""
struct ProjectedDensityOfStates{T} <: Vector{DensityOfStates{T}}
    states::Array{StateDensity{T},3}
    labels::Vector{InlineString15}
    ecut::T
    fermi::T
    hastotal::Bool
end

Base.has_offset_axes(p::ProjectedDensityOfStates) = !p.hastotal

"""
    has_total_dos(p::ProjectedDensityOfStates) -> Bool

Returns `true` if the total DOS data is present, in which case the `total` or 0 index contains the
corresponding `DensityOfStates` object.
"""
has_total_dos(p::ProjectedDensityOfStates) = p.hastotal

Base.size(p::ProjectedDensityOfStates) = reverse(size(p.states))

function Base.axes(p::ProjectedDensityOfStates)
    ax = reverse(axes(p.states))
    return Base.setindex(ax, ax[1] .- has_total_dos(p), 1)
end

function Base.getindex(p::ProjectedDensityOfStates, i::Int)
    !p.hastotal && iszero(i) && error("This dataset does not contain total density of states data.")
    return DensityOfStates(p.states[:,:,i+p.hastotal])
end

function Base.getindex(p::ProjectedDensityOfStates, s::AbstractString)
    s == "total" && return p[:,:,0]
    return DensityOfStates(p.states[:,:,findfirst(s, p.labels)])
end
