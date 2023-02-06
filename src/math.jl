"""
    Electrum._is_linearly_independent(vecs::AbstractMatrix{<:Number}) -> Bool
    Electrum._is_linearly_independent(vecs::AbstractVector{<:Number}...) -> Bool

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
    Electrum.reinterpret_index(sz::NTuple{D,<:Integer}, inds::Tuple) -> Tuple
    Electrum.reinterpret_index(g, inds::Tuple) -> Tuple

Converts indices provided in a call to `getindex()` to valid array indices of the backing field that
contains the data. This is intended for use with data structures that use zero-based indexing so
that data at the zero index corresponds to data at the origin.
"""
function reinterpret_index(sz::NTuple{D,<:Integer}, inds::Tuple) where D
    return ntuple(Val{length(inds)}()) do i
        inds[i] isa Colon ? Colon() : mod.(inds[i], sz[i]) .+ 1
    end
end

reinterpret_index(g, inds::Tuple) = reinterpret_index(size(g), inds)

"""
    Electrum.convert_to_transform(M, [dimensions]) -> AbstractMatrix{Int}

Converts a scalar, vector, or matrix that is intended to be used as a lattice transformation into
a transformation matrix.
"""
convert_to_transform(M::AbstractMatrix) = Int.(M)
convert_to_transform(v::AbstractVector) = diagm(Int.(v))
# Drop provided dimensionality for matrices/vectors, since it's already known from the input
convert_to_transform(x::AbstractVecOrMat, dimensions) = convert_to_transform(x)
# When dimensions are provided as integers, return a Matrix
# When provided as Val types, return an SMatrix, since we can dispatch on dimension
convert_to_transform(n::Real, dimensions::Integer) = diagm(fill(Int(n), dimensions))
convert_to_transform(n::Real, ::Val{D}) where D = diagm(fill(Int(n), SVector{D}))
convert_to_transform(U::UniformScaling, dimensions::Integer) = Int.(U(dimensions))
convert_to_transform(U::UniformScaling, ::Val{D}) where D = SMatrix{D,D,Int}(U)
