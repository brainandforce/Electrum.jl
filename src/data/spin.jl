"""
    SpinRange{M} <: AbstractUnitRange{Rational{Int}}

Represents a valid range of spin states corresponding with multiplicity `M`. These are equivalent to 
`UnitRange{Rational{Int}}` objects of the form `-(M - 1)//2:(M - 1)//2`.

These types are singleton types with a numeric `Int` parameter, allowing for its use as a dispatch
type.

# Examples
```
julia> SpinRange{3}() == -1:1
true

julia> SpinRange{4}() == -3//2:3//2
true
```
"""
struct SpinRange{M} <: AbstractUnitRange{Rational{Int}}
    SpinRange{M}() where M = (@assert M > 0 "M must be a positive integer."; new{Int(M)}())
end

SpinRange(M) = SpinRange{Int(M)}()

Base.size(::SpinRange{M}) where M = M
Base.first(::SpinRange{M}) where M = -(M-1)//2
Base.last(::SpinRange{M}) where M = (M-1)//2
Base.getindex(s::SpinRange{M}, i::Integer) where M = first(s) + (i - 1)

Base.UnitRange(::SpinRange) = first(s):last(s)
Base.convert(T::Type{<:AbstractVector}, s::SpinRange) = convert(T, UnitRange(s))
Base.show(io::IO, s::SpinRange) = print(io, typeof(s), "()")

function Base.show(io::IO, ::MIME"text/plain", s::SpinRange)
    show(io, s)
    print(io, " (equivalent to ", UnitRange(s), ")")
end
