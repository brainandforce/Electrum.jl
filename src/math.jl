"""
    Electrum._is_linearly_independent(M::AbstractMatrix) -> Bool
    Electrum._is_linearly_independent(vecs::AbstractVector...) -> Bool

Determines whether a set of vectors is linearly independent.

This function always returns `false` if the first dimension of the matrix is less than the second.
"""
is_linearly_independent(M::AbstractMatrix) = !isless(size(M)...) && rank(M) == minimum(size(M))
is_linearly_independent(vecs::AbstractVector...) = is_linearly_independent(hcat(vecs...))

"""
    Electrum.reinterpret_index(sz::NTuple{D,<:Integer}, i) -> typeof(i)
    Electrum.reinterpret_index(g, i) -> typeof(i)

Converts indices provided in a call to `getindex()` to valid array indices of the backing field that
contains the data. This is intended for use with data structures that use zero-based indexing so
that data at the zero index corresponds to data at the origin.
"""
function reinterpret_index(sz::NTuple{D,<:Integer}, inds::Tuple) where D
    return ntuple(Val{length(inds)}()) do i
        inds[i] isa Colon ? Colon() : mod.(inds[i], sz[i]) .+ 1
    end
end

function reinterpret_index(sz::NTuple{D,<:Integer}, i::CartesianIndex{D}) where D
    return CartesianIndex(mod.(i.I, sz) .+ 1)
end

reinterpret_index(g, i) = reinterpret_index(size(g), i)

"""
    Electrum.convert_to_transform(M, [dimensions]) -> AbstractMatrix{Int}

Converts a scalar, vector, or matrix that is intended to be used as a lattice transformation into
a transformation matrix.
"""
convert_to_transform(M::AbstractMatrix) = Int.(M)
convert_to_transform(v::AbstractVector) = diagm(Int.(v))
convert_to_transform(v::AbstractVector, ::Val{D}) where D = diagm(SVector{D,Int}(v))
# Drop provided dimensionality for matrices/vectors, since it's already known from the input
convert_to_transform(x::AbstractVecOrMat, dimensions) = convert_to_transform(x)
# When dimensions are provided as integers, return a Matrix
# When provided as Val types, return an SMatrix, since we can dispatch on dimension
convert_to_transform(n::Real, dimensions::Integer) = diagm(fill(Int(n), dimensions))
convert_to_transform(n::Real, ::Val{D}) where D = diagm(fill(Int(n), SVector{D}))
convert_to_transform(U::UniformScaling, dimensions::Integer) = Int.(U(dimensions))
convert_to_transform(U::UniformScaling, ::Val{D}) where D = SMatrix{D,D,Int}(U)

#---FFT ranges-------------------------------------------------------------------------------------#
"""
    FFTBins{D} <: AbstractArray{CartesianIndex{D},D}

An iterable object defining a range of integer FFT bins. This can be used to convert a Cartesian
array index to an FFT bin index, or vice versa.

The outputs use the convention where frequencies at or above the Nyquist frequency for that
dimension are negative, matching the output of `FFTW.fftfreq`.

```jldoctest
julia> FFTBins(4)
4-element FFTIndices{1}:
 CartesianIndex(0,)
 CartesianIndex(1,)
 CartesianIndex(-2,)
 CartesianIndex(-1,)

julia> FFTBins(3, 3)
3×3 FFTIndices{2}:
 CartesianIndex(0, 0)   CartesianIndex(0, 1)   CartesianIndex(0, -1)
 CartesianIndex(1, 0)   CartesianIndex(1, 1)   CartesianIndex(1, -1)
 CartesianIndex(-1, 0)  CartesianIndex(-1, 1)  CartesianIndex(-1, -1)
```
"""
struct FFTBins{D} <: AbstractArray{CartesianIndex{D},D}
    size::NTuple{D,Int}
    FFTBins(x::NTuple{D}) where D = new{D}(Int.(x))
end

FFTBins(x::Number...) = FFTBins(x)
FFTBins(a::AbstractArray) = FFTBins(size(a))

Base.axes(r::FFTBins) = Base.OneTo.(r.size)
Base.size(r::FFTBins) = r.size
Base.IndexStyle(::Type{<:FFTBins}) = IndexLinear()

function Base.getindex(r::FFTBins{D}, i::CartesianIndex{D}) where D
    return CartesianIndex(mod.(Tuple(i) .+ div.(size(r), 2) .- 1, size(r)) .- div.(size(r), 2))
end

Base.getindex(r::FFTBins, i::Integer...) = r[CartesianIndex(i)]
Base.getindex(r::FFTBins, i::Integer) = r[CartesianIndices(r)[i]]

Base.iterate(r::FFTBins, i = 1) = i in eachindex(r) ? (r[i], i+1) : nothing

"""
    FFTLength <: AbstractVector{Int}

The one-dimensional counterpart to `FFTBins`, supporting only a single dimension. Its elements are
plain `Int` types rather than the `CartesianIndex` of `FFTBins`.

In essence, it serves as a counterpart to `Base.OneTo` for FFT bins.

```jldoctest
julia> Electrum.FFTLength(4)
4-element FFTLength:
 0
 1
 -2
 -1
```
"""
struct FFTLength <: AbstractVector{Int}
    size::Int
    FFTLength(x::Number) = new(Int(x))
end

FFTLength(v::AbstractVector) = FFTLength(length(v))

Base.axes(r::FFTLength) = (Base.OneTo(r.size),)
Base.size(r::FFTLength) = (r.size,)
Base.IndexStyle(::Type{<:FFTLength}) = IndexLinear()

function Base.getindex(r::FFTLength, i::Integer)
    return mod(i + div(length(r), 2) - 1, length(r)) - div(length(r), 2)
end

Base.getindex(r::FFTLength, i::CartesianIndex{1}) = r[only(i.I)]

Base.iterate(r::FFTLength, i = 1) = i in eachindex(r) ? (r[i], i+1) : nothing

# For Julia 1.6 compatibility: must use keyword arguments
Base.sort(r::FFTLength) = range(minimum(r), stop = maximum(r))
