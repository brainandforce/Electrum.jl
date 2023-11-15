"""
    AbstractEnergyData{T<:Real}

Supertype for all data associated with energy values, with all quantities assumed to be in Hartree
atomic units and of numeric type `T`.
"""
abstract type AbstractEnergyData{T<:Real}
end

Base.eltype(::Type{AbstractEnergyData{T}}) where T = T

"""
    energy(e::AbstractEnergyData{T}) -> T

Returns the energy in hartree.

Types descending from `AbstractEnergyData{T}` with a field named `energy` will automatically have
this function defined.
"""
energy(e::AbstractEnergyData) = e.energy

"""
    occupancy(e::AbstractEnergyData{T}) -> T

Returns occupancy data associated with energy data. Usually, occupancy values range from 0 to 2 in
restricted calculations (no spin) or 0 to 1 in unrestricted calcuations.

Types descending from `AbstractEnergyData{T}` with a field named `occupancy` will automatically have
this function defined. For those which do not store occupancy data, the error message `"No occupancy
data is associated with this type"` is automatically thrown.
"""
function occupancy(e::AbstractEnergyData)
    hasproperty(e, :occupancy) || error("No occupancy data is associated with this type.")
    return e.occupancy
end

"""
    EnergyOccupancy{T} <: AbstractEnergyData{T}

A data structure consisting of a pair of an energy value and an occupancy number, both of type `T`.

These types serve as the elements of `BandStructure{T}`.
"""
struct EnergyOccupancy{T} <: AbstractEnergyData{T}
    energy::T
    occupancy::T
end

"""
    StateDensity{T} <: AbstractEnergyData{T}

A set of energy data, occupancy values, and the density of states at that energy, all of type `T`.

These types serve as the elements of `DensityOfStates{T}` and `ProjectedDensityOfStates{T}`.
"""
struct StateDensity{T} <: AbstractEnergyData{T}
    energy::T
    occupancy::T
    density::T
end

(T::Type{<:EnergyOccupancy})(s::StateDensity) = T(energy(s), occupancy(s))
Base.convert(T::Type{<:EnergyOccupancy}, s::StateDensity) = T(s)

"""
    density(s::StateDensity{T}) -> T

Returns the density of states at the energy in the `StateDensity` object.
"""
density(s::StateDensity) = s.density

# TODO: deprecate this type alias in favor of BandStructure{T}
"""
    EnergiesOccupancies{T,N} <: AbstractArray{EnergyOccupancy{T},N}

Type alias for `Array{EnergyOccupancy{T},N}`. Data structures `S` which contain `EnergyOccupancy`
values in a collection should define the constructor `EnergiesOccupancies(::S)`.
"""
const EnergiesOccupancies{T,N} = Array{EnergyOccupancy{T},N}

"""
    energies(x) -> Array{<:Real}

Returns the energy data associated with a collection of `EnergyOccupancy{T}` objects. By default,
this falls back to `energy.(EnergyOccupancy(x))`.
"""
energies(a::AbstractArray{<:AbstractEnergyData}) = energy.(a)
energies(x) = energy.(EnergiesOccupancies(x))

"""
    occupancies(x) -> Array{<:Real}

Returns the occupancy data associated with a collection of `EnergyOccupancy{T}` objects. By default,
this falls back to `occupancy.(EnergyOccupancy(x))`.
"""
occupancies(a::AbstractArray{<:AbstractEnergyData}) = occupancy.(a)
occupancies(x) = occupancy.(EnergiesOccupancies(x))

"""
    densities(x) -> Array{<:Real}

Returns the state density data associated with a collection of `StateDensity{T}` objects.
"""
densities(a::AbstractArray{<:StateDensity}) = density.(a)

"""
    min_energy(x) -> Real

Returns the minimum energy in a collection of EnergyOccupancy data or an object containing such
data.
"""
min_energy(a::AbstractArray{<:AbstractEnergyData}) = minimum(energy(eo) for eo in a)
min_energy(x) = min_energy(EnergiesOccupancies(x))

"""
    max_energy(x) -> Real

Returns the maximum energy in a collection of EnergyOccupancy data or an object containing such 
data.
"""
max_energy(a::AbstractArray{<:AbstractEnergyData}) = maximum(energy(eo) for eo in a)
max_energy(x) = max_energy(EnergiesOccupancies(x))

"""
    min_occupancy(x) -> Real

Returns the minimum occupancy in a collection of EnergyOccupancy data or an object containing such
data. For most raw calculation data, this should return zero if unoccupied states were considered.
"""
min_occupancy(a::AbstractArray{<:AbstractEnergyData}) = minimum(occupancy(eo) for eo in a)
min_occupancy(x) = min_occupancy(EnergiesOccupancies(x))

"""
    max_occupancy(x) -> Real

Returns the maximum occupancy in a collection of EnergyOccupancy data or an object containing such
data. For a restricted calculation (no explicit treatment of spin), this is usually around 2, and
for an unrestricted calculation (explicit spin treatment) this is usually around 1.

In many cases, you may want to determine the maximum possible occupancy value, not the maximum in
the dataset, in which case, you should use `round(Int, max_occupancy(a), RoundUp)`.
"""
max_occupancy(a::AbstractArray{<:AbstractEnergyData}) = maximum(occupancy(eo) for eo in a)
max_occupancy(x) = max_occupancy(EnergiesOccupancies(x))

"""
    fermi(x) -> Real

Estimates the Fermi energy using the provided energy and occupancy data by interpolating the data
between the occupancies nearest half of the maximum occupancy.
"""
function fermi(a::AbstractArray{<:AbstractEnergyData})
    # TODO: use Levenberg-Marquardt fitting
    maxocc = round(Int, max_occupancy(a))
    @assert maxocc in 1:2 "The calculated maximum occupancy was $maxocc."
    sorted = sort!(vec(a), by=energy)
    # Find the index of the last occupancy that's greater than half of maxocc
    ind = findlast(x -> occupancy(x) > maxocc/2, sorted)
    # Perform a linear interpolation between O[ind] and O[ind+1] 
    approx = (occupancy(sorted[ind]) - maxocc/2) / mapreduce(occupancy, -, sorted[ind:ind+1])
    return (energy(sorted[ind]) * approx) + (energy(sorted[ind+1]) * (1 - approx))
end

fermi(x) = fermi(EnergiesOccupancies(x))
