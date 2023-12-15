"""
    Electrum.ZeroTo{T<:Integer} <: AbstractUnitRange{T}

A range that is guaranteed to begin with `zero(T)`. This is similar to `Base.OneTo`, but intended 
for use with arrays that index from zero.
"""
struct ZeroTo{T<:Integer} <: AbstractUnitRange{T}
    stop::T
    function ZeroTo{T}(stop::Integer) where T
        if stop < 0 && !(T <: Signed)
            throw(ArgumentError("Final value must be positive for element type $T."))
        end
        return new(max(zero(T) - oneunit(T), stop))
    end
end

function ZeroTo{T}(r::AbstractRange) where T
    iszero(first(r)) || throw(ArgumentError("First element must be 0 (got $(first(r)))."))
    isone(step(r)) || throw(ArgumentError("Step size must be 1 (got $(step(r)))."))
    return ZeroTo{T}(last(r))
end

ZeroTo(i::Integer) = ZeroTo{typeof(i)}(i)
ZeroTo(r::AbstractRange) = ZeroTo{eltype(r)}(r)

Base.size(z::ZeroTo) = tuple(z.stop + oneunit(eltype(z)))
Base.first(z::ZeroTo) = zero(eltype(z))

Base.@propagate_inbounds Base.getindex(z::ZeroTo, i::Int) = (@boundscheck checkbounds(z, i); i - 1)

Base.UnitRange(z::ZeroTo) = first(z):last(z)
Base.show(io::IO, z::ZeroTo{Int}) = print(io, ZeroTo, '(', z.stop, ')')
Base.show(io::IO, z::ZeroTo) = print(io, typeof(z), '(', z.stop, ')')

#---New implementation of data grids associated with lattices--------------------------------------#
"""
    LatticeData{D,T,M<:LatticeDataMap{<:BySpace,D},A<:AbstractArray{T,D}} <: AbstractArray{T,D}

Maps an array of type `A` onto regularly spaced points of a lattice with a mapping of type `M`.
The array and lattice must share the same dimension parameter `D`.

# Indexing

Arrays wrapped by this data structure use periodic indexing, setting the first Cartesian index to
`zero(CartesianIndex{D})`, and first linear index to 0.

The `IndexStyle` matches that of `A`.
"""
struct LatticeData{D,T,M<:LatticeDataMap{<:BySpace,D},A<:AbstractArray{T,D}} <: AbstractArray{T,D}
    data::A
    map::M
end

"""
    RealLatticeData{D,T,X,A} (alias for LatticeData{D,T,LatticeDataMap{ByRealSpace,D,X},A})

An array `A` mapped onto regularly spaced points associated with a real space lattice.

For more information, see [LatticeData](@ref).
"""
const RealLatticeData{D,T,X,A} = LatticeData{D,T,LatticeDataMap{ByRealSpace,D,X},A}

"""
    ReciprocalLatticeData{D,T,X,A} (alias for LatticeData{D,T,LatticeDataMap{ByRealSpace,D,X},A})

An array `A` mapped onto regularly spaced points associated with a reciprocal space lattice.

For more information, see [LatticeData](@ref).
"""
const ReciprocalLatticeData{D,T,X,A} = LatticeData{D,T,LatticeDataMap{ByReciprocalSpace,D,X},A}

#=
const LatticeArray{M,D,T} = LatticeData{D,T,M,Array{T,D}}
const RealLatticeArray{X,D,T} = LatticeArray{LatticeDataMap{ByRealSpace,D,X},D,T}
const ReciprocalLatticeArray{X,D,T} = LatticeArray{LatticeDataMap{ByRealSpace,D,X},D,T}
=#

basis(l::LatticeData) = basis(l.map)
BySpace(::Type{<:LatticeData{D,M}}) where {D,M} = BySpace(M)
LatticeDataMap(l::LatticeData) = l.map

Base.size(l::LatticeData) = size(l.data)
Base.axes(l::LatticeData) = ZeroTo.(size(l) .- 1)

# Index style matches that of the wrapped array
Base.IndexStyle(::Type{T}) where T<:LatticeData = IndexStyle(fieldtype(T, :data))

Base.to_indices(l::LatticeData, I::Tuple) = map((x,i) -> mod.(x,i), to_indices(l.data, I), size(l))

Base.getindex(l::LatticeData, i...) = @inbounds getindex(l.data, to_indices(l, i)...)
Base.setindex!(l::LatticeData, x, i...) = @inbounds setindex!(l.data, x, to_indices(l, i)...)

#---Broadcasting-----------------------------------------------------------------------------------#
"""
    Electrum.is_broadcast_compatible(l::LatticeData...)

Checks that two lattices are compatible for a broadcasted operation by 
"""
is_broadcast_compatible(a::LatticeData, b::LatticeData) = LatticeDataMap(a) == LatticeDataMap(b)
is_broadcast_compatible(l::LatticeData...) = reduce(is_broadcast_compatible, l)

#---Mathematical operations------------------------------------------------------------------------#
