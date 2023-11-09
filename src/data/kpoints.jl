"""
    KPoint{D} <: DenseVector{Float64}

Stores a k-point with an associated weight that corresponds to the number of symmetry-equivalent
k-points, stored as an integer.
"""
struct KPoint{D} <: DenseVector{Float64}
    point::SVector{D,Float64}
    weight::Int
    KPoint(pt::StaticVector{D,<:Real}, wt::Integer = 1) where D = new{D}(pt, wt)
end

KPoint(x::Real...; weight::Integer = 1) = KPoint(SVector(x), weight)
KPoint{D}(pt::AbstractVector{<:Real}, wt::Integer = 1) where D = KPoint(SVector{D,Float64}(pt), wt)

Base.hash(k::KPoint, h::UInt) = hash(k.point, hash(k.weight, h))
Base.:(==)(k1::KPoint, k2::KPoint) = k1.point == k2.point && k1.weight == k2.weight

Base.length(::Type{KPoint{D}}) where D = D
Base.size(::Type{KPoint{D}}) where D = (D,)
Base.size(k::KPoint) = size(k.point)
Base.axes(::Type{KPoint{D}}) where D = (Base.OneTo(D),)
Base.axes(k::KPoint) = axes(k.point)

Base.IndexStyle(::Type{<:KPoint}) = IndexLinear()
Base.getindex(k::KPoint, i) = k.point[i]

Base.convert(T::Type{<:StaticVector}, k::KPoint) = convert(T, k.point)::T
Base.convert(T::Type{<:Vector}, k::KPoint) = convert(T, k.point)::T
Base.convert(T::Type{<:KPoint}, v::AbstractVector{<:Real}) = T(v, 1)
Base.convert(T::Type{<:KPoint}, k::KPoint) = k::T

Base.zero(::Type{KPoint{D}}) where D = KPoint(zero(SVector{D,Float64}))

"""
    weight(k::KPoint) -> Int

Returns the weight associated with a k-point.
"""
weight(k::KPoint) = k.weight

# TODO: can we implement a remainder that excludes -0.5?
"""
    truncate(k::KPoint) -> KPoint

Moves a k-point so that its values lie within the range [-1/2, 1/2]. The weight is preserved.
"""
Base.truncate(k::KPoint) = KPoint(rem.(k.point, 1, RoundNearest), k.weight)

#---Generated lists of k-points--------------------------------------------------------------------#
"""
    KPointMesh{D} <: AbstractVector{KPoint{D}}

Contains a list of k-points associated with a matrix describing the mesh that was used to generate
the points, and its shift off the Γ point (origin). If the mesh used to generate the points is
unknown, it will be set to the zero matrix of dimension `D`.

A `KPointMesh` can be indexed as if it were an ordinary `Vector{KPoint{D}}`.
"""
struct KPointMesh{D} <: AbstractVector{KPoint{D}}
    points::Vector{KPoint{D}}
    grid::SMatrix{D,D,Int}
    shift::SVector{D,Float64}
    function KPointMesh(
        points::AbstractVector{KPoint{D}}, 
        grid::StaticMatrix{D,D,<:Integer} = zeros(SMatrix{D,D,Int}),
        shift::StaticVector{D,<:Real} = zeros(SVector{D,Float64})
    ) where D
        return new{D}(points, grid, shift)
    end
end

function KPointMesh(
    points::AbstractVector,
    grid::StaticMatrix{D,D,<:Integer} = zeros(SMatrix{D,D,Int}),
    shift::StaticVector{D,<:Real} = zeros(SVector{D,Float64})
) where D
    return KPointMesh(KPoint{D}.(points), grid, shift)
end

function KPointMesh{D}(
    points::AbstractVector,
    grid::AbstractMatrix{<:Integer} = zeros(SMatrix{D,D,Int}),
    shift::AbstractVector{<:Real} = zeros(SVector{D,Float64})
) where D
    return KPointMesh(points, SMatrix{D,D,Int}(grid), SVector{D,Float64}(shift))
end

function Base.:(==)(k1::KPointMesh, k2::KPointMesh)
    return k1.points == k2.points && k1.grid == k2.grid && k1.shift == k2.shift
end

Base.size(k::KPointMesh) = size(k.points)
Base.axes(k::KPointMesh) = axes(k.points)

Base.IndexStyle(::Type{<:KPointMesh}) = IndexLinear()
Base.getindex(k::KPointMesh, i) = k.points[i]
Base.setindex!(k::KPointMesh, x, i) = setindex!(k.points, x, i)

Base.convert(T::Type{Vector{<:KPoint}}, k::KPointMesh) = k.points::T
Base.convert(T::Type{<:KPointMesh}, v::AbstractVector{<:KPoint}) = KPointMesh(v)::T
Base.convert(T::Type{<:KPointMesh}, k::KPointMesh) = k::T

"""
    nkpt(k::KPointMesh) -> Int

Counts the number of k-points (equivalent to `length(k)`). This function can be defined for custom
types that contain a `KPointMesh`.
"""
nkpt(k::KPointMesh) = length(k)

"""
    nkpt(x) -> Int

Counts the number of explicitly enumerated k-points associated with an object containing a
`KPointMesh`.

By default, this function returns `length(KPointMesh(x))`, so defining `KPointMesh(x::T)` for a
type `T` containing a `KPointMesh` will automatically define `nkpt(x::T)`.
"""
nkpt(x) = length(KPointMesh(x))

#---Energy/occupancy pairs-------------------------------------------------------------------------#
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
energies(x) = energy.(EnergyOccupancy(x))

"""
    occupancies(x) -> Array{<:Real}

Returns the occupancy data associated with a collection of `EnergyOccupancy{T}` objects. By default,
this falls back to `occupancy.(EnergyOccupancy(x))`.
"""
occupancies(a::AbstractArray{<:EnergyOccupancy}) = occupancy.(a)
occupancies(x) = occupancy.(EnergyOccupancy(x))

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
