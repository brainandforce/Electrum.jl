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
