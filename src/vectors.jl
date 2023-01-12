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
    Xtal.reinterpret_index(sz::NTuple{D,<:Integer}, inds::Tuple) -> Tuple
    Xtal.reinterpret_index(g, inds::Tuple) -> Tuple

Converts indices provided in a call to `getindex()` to valid array indices of the backing field
that contains the data. This is intended for use with data structures that use zero-based indexing
so that data at the zero index corresponds to data at the origin.
"""
function reinterpret_index(sz::NTuple{D,<:Integer}, inds::Tuple) where D
    return ntuple(Val{length(inds)}()) do i
        inds[i] isa Colon ? Colon() : mod.(inds[i], sz[i]) .+ 1
    end
end

reinterpret_index(g, inds::Tuple) = reinterpret_index(size(g), inds)
