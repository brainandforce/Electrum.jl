"""
    Xtal._is_linearly_independent(vecs::AbstractMatrix{<:Number}) -> Bool
    Xtal._is_linearly_independent(vecs::AbstractVector{<:Number}...) -> Bool

Determines whether a set of vectors is linearly independent.

This function always returns `false` if the first dimension of the matrix is less than the second.
"""
function is_linearly_independent(vecs::AbstractMatrix{<:Number})
    size(vecs)[1] < size(vecs)[2] || return false
    return LinearAlgebra.rank(vecs) == minimum(size(vecs))
end

# Same thing but for sets of vectors
function is_linearly_independent(vecs::AbstractVector{<:Number}...)
    return is_linearly_independent(hcat(vecs...))
end

"""
    Xtal.fbin(i::Integer) -> Vector{Int}
    Xtal.fbin(v::AbstractVector) -> Vector{Int}

Returns the coefficients corresponding to the smallest positive and negative frequencies 
associated with each frequency bin in a discrete Fourier transform of a vector. For 
multidimensional arrays, use `fbins()`.

For vectors of even length, the maximum frequency element can be interpreted as either a positive
or negative value; it is selected to be positive.

# Examples

```jldoctest
julia> fbin(6)
6-element Vector{Int64}:
  0
  1
  2
  3
 -2
 -1
```
"""
fbin(i::Integer) = vcat(collect(0:div(i, 2)), collect(-div(i-1, 2):-1))
fbin(v::AbstractVector) = fbin(length(v))
fbin(::SVector{N,<:Any}) where N = SVector{N,Int}(fbin(N))

"""
    Xtal.fbins(i::Vararg{T,N}) -> Array{NTuple{N,Int},N}
    Xtal.fbins(x::AbstractArray{T,N}) -> Array{NTuple{N,Int},N}

Returns the coefficients corresponding to the smallest positive and negative frequencies 
associated with each frequency bin in a discrete Fourier transform. This is intended to be used
with multidimensional arrays, where there are multiple coefficients that need to be tracked. 
Vectors can use `fbin()` to get a `Vector{Int}` that should be easier to work with.

For vectors of even length, the maximum frequency element can be interpreted as either a positive
or negative value; it is selected to be positive.

# Examples

```jldoctest
julia> fbins(6,5)
6×5 Matrix{Tuple{Int64, Int64}}:
 (0, 0)   (0, 1)   (0, 2)   (0, -2)   (0, -1)
 (1, 0)   (1, 1)   (1, 2)   (1, -2)   (1, -1)
 (2, 0)   (2, 1)   (2, 2)   (2, -2)   (2, -1)
 (3, 0)   (3, 1)   (3, 2)   (3, -2)   (3, -1)
 (-2, 0)  (-2, 1)  (-2, 2)  (-2, -2)  (-2, -1)
 (-1, 0)  (-1, 1)  (-1, 2)  (-1, -2)  (-1, -1)
```   
"""
fbins(i::Integer...) = collect(Iterators.product((fbin(x) for x in i)...))
fbins(a::AbstractArray) = fbins(size(a)...)
fbins(s::SArray{S,<:Any,N}) where {S,N} = SArray{S,NTuple{N,Int}}(fbins(size(s)...))

#=
"""
    _to_vectors(a::AbstractArray{T,N}) where {T,N} -> Vector{...Vector{T}}

Converts an array to a nest of vectors within vectors.
"""
function _to_vectors(a::AbstractArray{T,N}) where {T,N}
    
end
=#