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

"""
    nspin(b::BandStructure) -> Int

Returns the number of spins associated with `b` (usually 1 or 2).
"""
nspin(b::BandStructure) = size(b, 1)

"""
    nkpt(b::BandStructure) -> Int

Returns the number of k-points associated with `b`.
"""
function nkpt(b::BandStructure)
    @assert size(b, 2) == length(b.kpoints) "Mismatch in expected number of k-points"
    return size(b, 2)
end

"""
    nband(b::BandStructure) -> Int

Returns the number of bands associated with `b`.
"""
nband(b::BandStructure) = size(b, 3)

KPointMesh(b::BandStructure) = b.kpoints

"""
    ecut(b::BandStructure{D,T}) -> Union{T,Missing}

Returns the energy cutoff of `b`. By default, this returns `b.ecut`, unless its value is `NaN`, in
which case it returns `missing`.
"""
ecut(b::BandStructure) = ifelse(isnan(b.ecut), missing, b.ecut)

"""
    fermi(b::BandStructure{D,T}) -> Union{T,Missing}

Returns the Fermi energy of `b`. By default, this returns `b.fermi`, unless its value is `NaN`, in
which case it returns `missing`.
"""
fermi(b::BandStructure) = ifelse(isnan(b.fermi), missing, b.fermi)

#---Density of states------------------------------------------------------------------------------#

