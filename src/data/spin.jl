"""
    Multiplicity{M} <: AbstractUnitRange{Rational{Int}}

Represents a valid range of spin states corresponding with multiplicity `M`. These index as if they 
are `UnitRange{Rational{Int}}` objects of the form `-(M - 1)//2:(M - 1)//2`, and can be converted to
`UnitRange` objects with the `UnitRange` constructor.

These types are singleton types with a numeric `Int` parameter, allowing for its use as a type
parameter in other types which require information on the spin multiplicity or allowed spin states.

# Examples
```
julia> UnitRange(SpinRange{3}())
1:1

julia> SpinRange{4}() == -3//2:3//2
true
```
"""
struct Multiplicity{M} <: AbstractUnitRange{Rational{Int}}
    Multiplicity{M}() where M = (@assert M > 0 "M must be a positive integer."; new{Int(M)}())
end

"""
    Multiplicity(M)

Constructs `Multiplicity{M}()`. This function is not type-stable unless the value of `M` is known at
compile time, similar to `Val(M)`.
"""
Multiplicity(M) = Multiplicity{Int(M)}()

Base.size(::Multiplicity{M}) where M = M
Base.first(::Multiplicity{M}) where M = -(M-1)//2
Base.last(::Multiplicity{M}) where M = (M-1)//2

@inbounds function Base.getindex(s::Multiplicity, i::Integer)
    return i in eachindex(s) ? first(s) + (i - 1) : throw(BoundsError(s, i))
end

Base.UnitRange(s::Multiplicity) = first(s):last(s)
Base.show(io::IO, s::Multiplicity) = print(io, typeof(s), "()")
