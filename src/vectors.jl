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
