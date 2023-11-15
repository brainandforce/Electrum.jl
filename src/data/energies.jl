"""
    EnergyOccupancy{T<:Real}

A data structure consisting of a pair of an energy value and an occupancy number, both of type `T`.
Energies are assumed to be in Hartree. Occupancy values are not constrained, but will generally
range from 0 to 2 for the results of restricted calculations (no separate treatment of spins), or
from 0 to 1 for unrestricted calculations (wavefunctions with separate spins).
"""
struct EnergyOccupancy{T<:Real}
    energy::T
    occupancy::T
end

"""
    EnergiesOccupancies{T,N} <: AbstractArray{EnergyOccupancy{T},N}

Type alias for `Array{EnergyOccupancy{T},N}`. Data structures `S` which contain `EnergyOccupancy`
values in a collection should define the constructor `EnergiesOccupancies(::S)`.
"""
const EnergiesOccupancies{T,N} = Array{EnergyOccupancy{T},N}

"""
    energy(eo::EnergyOccupancy{T}) -> T

Returns the energy value in an `EnergyOccupancy`.
"""
energy(eo::EnergyOccupancy) = eo.energy

"""
    occupancy(eo::EnergyOccupancy{T}) -> T

Returns the occupancy value in an `EnergyOccupancy`.
"""
occupancy(eo::EnergyOccupancy) = eo.occupancy

"""
    energies(x) -> Array{<:Real}

Returns the energy data associated with a collection of `EnergyOccupancy{T}` objects. By default,
this falls back to `energy.(EnergyOccupancy(x))`.
"""
energies(a::AbstractArray{<:EnergyOccupancy}) = energy.(a)
energies(x) = energy.(EnergiesOccupancies(x))

"""
    occupancies(x) -> Array{<:Real}

Returns the occupancy data associated with a collection of `EnergyOccupancy{T}` objects. By default,
this falls back to `occupancy.(EnergyOccupancy(x))`.
"""
occupancies(a::AbstractArray{<:EnergyOccupancy}) = occupancy.(a)
occupancies(x) = occupancy.(EnergiesOccupancies(x))

"""
    min_energy(x) -> Real

Returns the minimum energy in a collection of EnergyOccupancy data or an object containing such
data.
"""
min_energy(a::AbstractArray{<:EnergyOccupancy}) = minimum(energy(eo) for eo in a)
min_energy(x) = min_energy(EnergiesOccupancies(x))

"""
    max_energy(x) -> Real

Returns the maximum energy in a collection of EnergyOccupancy data or an object containing such 
data.
"""
max_energy(a::AbstractArray{<:EnergyOccupancy}) = maximum(energy(eo) for eo in a)
max_energy(x) = max_energy(EnergiesOccupancies(x))

"""
    min_occupancy(x) -> Real

Returns the minimum occupancy in a collection of EnergyOccupancy data or an object containing such
data. For most raw calculation data, this should return zero if unoccupied states were considered.
"""
min_occupancy(a::AbstractArray{<:EnergyOccupancy}) = minimum(occupancy(eo) for eo in a)
min_occupancy(x) = min_occupancy(EnergiesOccupancies(x))

"""
    max_occupancy(x) -> Real

Returns the maximum occupancy in a collection of EnergyOccupancy data or an object containing such
data. For a restricted calculation (no explicit treatment of spin), this is usually around 2, and
for an unrestricted calculation (explicit spin treatment) this is usually around 1.

In many cases, you may want to determine the maximum possible occupancy value, not the maximum in
the dataset, in which case, you should use `round(Int, max_occupancy(a), RoundUp)`.
"""
max_occupancy(a::AbstractArray{<:EnergyOccupancy}) = maximum(occupancy(eo) for eo in a)
max_occupancy(x) = max_occupancy(EnergiesOccupancies(x))

"""
    fermi(x) -> Real

Estimates the Fermi energy using the provided energy and occupancy data by interpolating the data
between the occupancies nearest half of the maximum occupancy.
"""
function fermi(a::AbstractArray{<:EnergyOccupancy})
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
