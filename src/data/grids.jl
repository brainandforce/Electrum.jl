"""
    Electrum.correct_shift(v::AbstractVector{<:Real})

Ensures that the components of a shift vector used in the construction of a `DataGrid` remain close
to zero - in other words, that they do not exceed 0.5.
"""
correct_shift(v::AbstractVector{<:Real}) = v - round.(v)
# This check is already done for KPoint
correct_shift(k::KPoint) = k

#---Struct definition------------------------------------------------------------------------------#
"""
    Electrum.DataGrid{D,B<:LatticeBasis,S<:AbstractVector{<:Real},T} <: AbstractArray{T,D}

Stores a grid of values of type `T` defined in a `D`-dimensional crystal lattice basis of type `B`
with a shift parameter of type `S`.

By convention, indexing of `DataGrid` is zero-based. This convention is used so that the first entry
corresponds to data at the origin can be indexed with zeros. However, `getindex()` is implemented
such that the dataset may be indexed by any integer, with modulo math used to convert to an index
within the grid. This is done with `Electrum.reinterpret_index()`.

Linear indexing is one-based like that of the underlying array.

For convenience, the aliases `RealDataGrid` and `ReciprocalDataGrid` are provided and are defined
below:

    const RealDataGrid{D} = DataGrid{RealBasis{D,Float64},SVector{D,Float64},D}
    const ReciprocalDataGrid{D} = DataGrid{ReciprocalBasis{D,Float64},KPoint{D},D}

Note that `ReciprocalDataGrid` uses a `KPoint{D}` to represent the shift, as opposed to an 
`SVector{D}`. This allows for the concurrent storage of a weight along with the shift, which may be
relevant for wavefunctions which exploit the symmetry of the k-point mesh.
"""
struct DataGrid{D,B<:LatticeBasis,S<:AbstractVector{<:Real},T} <: AbstractArray{T,D}
    data::Array{T,D}
    basis::B
    shift::S
    function DataGrid{D,B,S,T}(
        data::AbstractArray,
        basis::AbstractMatrix{<:Real},
        shift::AbstractVector{<:Real} = zero(S)
    ) where {D,B,S,T}
        return new(data, basis, correct_shift(shift))
    end
end

Base.has_offset_axes(g::DataGrid) = true

Base.show(io::IO, g::DataGrid) = print(io, typeof(g), (g.data, g.basis, g.shift))

const RealDataGrid{D,T} = DataGrid{D,RealBasis{D,Float64},SVector{D,Float64},T}
const ReciprocalDataGrid{D,T} = DataGrid{D,ReciprocalBasis{D,Float64},KPoint{D},T}

#---Constructors-----------------------------------------------------------------------------------#
"""
    Electrum.zero_shift(::Type{<:LatticeBasis{S,D}})

Generates a default zero shift vector whose type depends on the type of the basis vector.

This function is needed because `RealDataGrid` uses `SVector{D,Float64}` as its shift vector, and
`ReciprocalDataGrid` uses `KPoint{D}`. Depending on the basis vector type used in the constructor,
the shift vector type must be automatically determined if not supplied.
"""
zero_shift(::Type{<:RealBasis{D}}) where D = zero(SVector{D,Float64})
zero_shift(::Type{<:ReciprocalBasis{D}}) where D = zero(KPoint{D})
zero_shift(b::LatticeBasis) = zero_shift(typeof(b))

function DataGrid{D,B,S}(
    data::AbstractArray{T},
    basis::AbstractMatrix{<:Real},
    shift::AbstractVector{<:Real} = zero_shift(B)
) where {D,B,S,T}
    return DataGrid{D,B,S,T}(data, basis, shift)
end

function DataGrid{D,B}(
    data::AbstractArray{T},
    basis::AbstractMatrix{<:Real},
    shift::S = zero_shift(B)
) where {D,B,S,T}
    return DataGrid{D,B,S,T}(data, basis, shift)
end

function DataGrid{D}(data::AbstractArray{T}, basis::B, shift::S = zero_shift(B)) where {D,B,S,T}
    return DataGrid{D,B,S,T}(data, basis, shift)
end

function DataGrid(data::AbstractArray{T,D}, basis::B, shift::S = zero_shift(B)) where {D,B,S,T}
    return DataGrid{D,B,S,T}(data, basis, shift)
end

function RealDataGrid(
    data::AbstractArray{T,D},
    basis::AbstractMatrix{<:Real},
    shift::AbstractVector{<:Real} = zero(SVector{D,Float64})
) where {D,T}
    return DataGrid{D,RealBasis{D,Float64},SVector{D,Float64},T}(data, basis, shift)
end

function ReciprocalDataGrid(
    data::AbstractArray{T,D},
    basis::AbstractMatrix{<:Real},
    shift::AbstractVector{<:Real} = zero(KPoint{D})
) where {D,T}
    return DataGrid{D,ReciprocalBasis{D,Float64},KPoint{D},T}(data, basis, shift)
end

#---Initialize array with zeros--------------------------------------------------------------------#

function Base.zeros(
    ::Type{DataGrid{D,B,S,T}},
    basis::LatticeBasis,
    shift::AbstractVector{<:Real},
    dimensions::Vararg{Integer,D}
) where {D,B,S,T}
    return T(zeros(T, dimensions), basis, shift)
end

function Base.zeros(
    ::Type{DataGrid{D,B,S,T}},
    basis::LatticeBasis,
    dimensions::Vararg{Integer,D}
) where {D,B,S,T}
    return DataGrid{D,B,S,T}(zeros(T, dimensions), basis, zero(S))
end

#---Equality and hashing---------------------------------------------------------------------------#

function Base.:(==)(g::DataGrid, h::DataGrid)
    return all(getfield(g, f) == getfield(h, f) for f in fieldnames(DataGrid))
end

Base.hash(g::DataGrid, h::UInt) = reduce(xor, hash(getfield(g, f), h) for f in fieldnames(DataGrid))

#---Size, axes, indexing---------------------------------------------------------------------------#

Base.size(g::DataGrid) = size(g.data)
# For Julia 1.6 compatibility: must use keyword arguments
Base.axes(g::DataGrid{D}) where D = NTuple{D}(range(0, stop = x - 1) for x in size(g))

# Linear indexing should be defined automatically
Base.getindex(g::DataGrid, i...) = getindex(g.data, reinterpret_index(g, i)...)
Base.setindex!(g::DataGrid, x, i...) = setindex!(g.data, x, reinterpret_index(g, i)...)

Base.getindex(g::DataGrid, i::CartesianIndex) = getindex(g, i.I...)
Base.setindex!(g::DataGrid, x, i::CartesianIndex) = setindex!(g, x, i.I...)

Base.CartesianIndices(g::DataGrid) = CartesianIndices(axes(g))

#---Get basis vectors and grid shift---------------------------------------------------------------#

(T::Type{<:LatticeBasis})(g::DataGrid) = convert(T, basis(g))
shift(g::DataGrid) = g.shift

#---Volumes associated with data grids-------------------------------------------------------------#
"""
    volume(g::DataGrid) -> eltype(basis(g))

Gets the crystal volume associated with a `DataGrid`. Units are assumed to be bohr³ for `DataGrid`
types with a `RealBasis`, and rad³*bohr⁻³ for types with a `ReciprocalBasis`.
"""
volume(g::DataGrid) = volume(basis(g))

"""
    voxelsize(g::DataGrid) -> eltype(basis(g))

Gets the size of a single real space voxel associated with a data grid. Units are assumed to be
bohr³.
"""
voxelsize(g::DataGrid) = volume(RealBasis(g)) / length(g)

#---Unary mathematical operations------------------------------------------------------------------#

Base.:-(g::DataGrid) = DataGrid(-g.data, g.basis, g.shift)

Base.abs(g::DataGrid) = DataGrid(abs.(g.data), g.basis, g.shift)
Base.abs2(g::DataGrid) = DataGrid(abs2.(g.data), g.basis, g.shift)
Base.conj(g::DataGrid) = DataGrid(conj.(g.data), g.basis, g.shift)
Base.angle(g::DataGrid) = DataGrid(angle.(g.data) + 2π * (imag.(g.data) .< 0), g.basis, g.shift)

#---Collecting shared grid properties--------------------------------------------------------------#
"""
    Electrum.get_shared_properties(f, args::Tuple, [exception::Exception = ErrorException("")])
        -> eltype(f.(args))

Promotes and compares the results of applying function or functor `f` to a set of arguments. If all
elements of the resulting tuple are equal, the only unique element is returned. If not, an exception
is thrown.
"""
function get_shared_properties(f, args::Tuple, ex::Exception = ErrorException(""))
    p = promote(f.(args)...)
    return all(==(first(p)), p) ? first(p) : throw(ex)
end

"""
    Electrum.get_basis_shift(grids::Tuple{Vararg{Datagrid}})

Gets the shared lattice basis vectors and shift associated with a set of datagrids.
"""
function get_basis_shift(args::Tuple{Vararg{DataGrid}})
    return (
        get_shared_properties(basis, args, LatticeMismatch("Unequal lattice basis vectors.")),
        get_shared_properties(shift, args, LatticeMismatch("Incommensurate lattice shifts."))
    )
end

#---Broadcasting-----------------------------------------------------------------------------------#

function Base.similar(g::DataGrid, ::Type{T}, sz::Tuple{Vararg{Int}}) where T
    return DataGrid(Array{T}(undef, sz), g.basis, g.shift)
end

function Base.similar(g::DataGrid, ::Type{T}, ax::Tuple{UnitRange,Vararg{UnitRange}}) where T
    return DataGrid(Array{T}(undef, length.(ax)), g.basis, g.shift)
end

function Base.similar(g::DataGrid, ::Type{T}, ax::FFTBins) where T
    return DataGrid(Array{T}(undef, size(ax)), g.basis, g.shift)
end

Base.similar(g::DataGrid, ::Type{T}) where T = similar(g, T, size(g))
Base.similar(g::DataGrid, sz::Tuple) = similar(g, eltype(g), sz)
Base.similar(g::DataGrid) = similar(g, eltype(g), size(g))

"""
    DataGridStyle{D,B<:LatticeBasis,S<:AbstractVector{<:Real}} <: Broadcast.AbstractArrayStyle{D}

The broadcast style for `DataGrid` objects.
"""
struct DataGridStyle{D,B<:LatticeBasis,S<:AbstractVector{<:Real}} <: Broadcast.AbstractArrayStyle{D}
end

DataGridStyle{D,B,S}(::Val{D}) where {D,B,S} = DataGridStyle{D,B,S}()
# Mismatching dimensions should fail
function DataGridStyle{D,B,S}(::Val) where {D,B,S}
    throw(DimensionMismatch("Cannot broadcast a DataGrid to arrays of varying dimensions"))
end

Base.BroadcastStyle(::Type{<:DataGrid{D,B,S}}) where {D,B,S} = DataGridStyle{D,B,S}()

# Allow broadcasting with other arrays, but don't preserve basis/shift - return Array
Base.BroadcastStyle(::DataGridStyle{D,B,S}, x::Broadcast.DefaultArrayStyle{D}) where {D,B,S} = x
Base.BroadcastStyle(::DataGridStyle{D,B,S}, x::Broadcast.AbstractArrayStyle{D}) where {D,B,S} = x

# Promote the element types of the basis and shift
function Base.BroadcastStyle(::DataGridStyle{D,B,S}, ::DataGridStyle{D,C,T}) where {D,B,S,C,T}
    return DataGridStyle{D,promote_type(B,C),promote_type(S,T)}()
end

function Base.similar(bc::Broadcast.Broadcasted{DataGridStyle{D,B,S}}, T::Type) where {D,B,S}
    grids = filter(x -> x isa DataGrid, bc.args)
    (b, s) = get_basis_shift(grids)
    sz = get_shared_properties(size, grids, DimensionMismatch("Unequal grid sizes."))
    return DataGrid(Array{T}(undef, sz), b, s)
end

#---More mathematical operations-------------------------------------------------------------------#

function Base.:≈(g::DataGrid, h::DataGrid; kw...)
    return all(getfield(g, f; kw...) ≈ getfield(h, f; kw...) for f in fieldnames(DataGrid))
end

Base.:*(s::Number, g::DataGrid) = DataGrid(s * g.data, g.basis, g.shift)
Base.:*(g::DataGrid, s::Number) = DataGrid(g.data * s, g.basis, g.shift)

Base.:/(s::Number, g::DataGrid) = DataGrid(s / g.data, g.basis, g.shift)
Base.:/(g::DataGrid, s::Number) = DataGrid(g.data / s, g.basis, g.shift)
# Right division is never used, but it's here for convenience
Base.:\(s::Number, g::DataGrid) = DataGrid(s \ g.data, g.basis, g.shift)
Base.:\(g::DataGrid, s::Number) = DataGrid(g.data \ s, g.basis, g.shift)

Base.:+(g::DataGrid, h::DataGrid) = DataGrid(g.data + h.data, get_basis_shift((g, h))...)
Base.:-(g::DataGrid, h::DataGrid) = DataGrid(g.data - h.data, get_basis_shift((g, h))...)

#---FFTs-------------------------------------------------------------------------------------------#
"""
    Electrum.fftvol(g::DataGrid)

Returns the volume needed to normalize the Fourier transform of `g`.
"""
fftvol(g::DataGrid{D,<:RealBasis}) where D = det(basis(g).matrix) / length(g)
fftvol(g::DataGrid{D,<:ReciprocalBasis}) where D = det(basis(g).matrix / 2π)

"""
    fft(g::RealDataGrid) -> ReciprocalDataGrid
    fft(g::ReciprocalDataGrid) -> RealDataGrid

Performs a fast Fourier transform on the data in a `DataGrid`.
"""
FFTW.fft(g::DataGrid) = DataGrid(fft(g.data) * fftvol(g), inv(basis(g)))

"""
    ifft(g::RealDataGrid) -> ReciprocalDataGrid
    ifft(g::ReciprocalDataGrid) -> RealDataGrid

Performs an inverse fast Fourier transform on the data in a `DataGrid`.

The inverse FFT is normalized so that `ifft(fft(g)) ≈ g` (to within floating point error).
"""
FFTW.ifft(g::DataGrid) = DataGrid(bfft(g.data) * fftvol(g), inv(basis(g)))

"""
    fftfreq(g::DataGrid{D}) -> Array{SVector{D,Float64},D}

Returns the Fourier transform frequency bins for an `DataGrid`.

For real space data, the frequency bins will be angular wavenumbers, matching the 2π factors that
are introduced when transforming between a `RealBasis` and a `ReciprocalBasis`. The convention used
by `FFTW.fftfreq()` is also used: frequency bins at or above the Nyquist frequency will be negative.

For reciprocal space data, the frequencies are binned with the assumption that the lattice vectors
are given in angular wavenumbers, and they represent real space coordinates. The Nyquist frequency
convention is *not* used, so all elements will have positive indices.
"""
function FFTW.fftfreq(g::DataGrid, ::ByRealSpace)
    return map(i -> SVector(2π .* Tuple(i) ./ size(g)), FFTBins(g))
end

function FFTW.fftfreq(g::DataGrid, ::ByReciprocalSpace)
    return map(i -> SVector(Tuple(i) .* size(g) ./ 2π), CartesianIndices(g))
end

FFTW.fftfreq(g::DataGrid) = fftfreq(g, DataSpace(g))

#---Operations specific to `RealDataGrid`----------------------------------------------------------#
"""
    integrate([f::Function = identity], g::RealDataGrid{D,T}) -> T

Applies the function `f` elementwise to a datagrid, then integrates the grid across all voxels,
accounting for the voxel volumes.
"""
integrate(f, g::RealDataGrid) = sum(f.(g.data)) * voxelsize(g)
integrate(g::RealDataGrid) = sum(g.data) * voxelsize(g)

"""
    Electrum.pdev_kernel(g::RealDataGrid{D}, dim::Integer) -> Array{T,D}

Calculates the partial derivative of `g` along the dimension indexed by `dim`, and return the result
in a plain Julia array. This is implemented with a fast Fourier transform.

A specialized method exists to calculate the partial derivative for real-valued data, which uses
`FFTW.rfft()` and `FFTW.brfft()` to cut down on the size of the final transform.
"""
function pdev_kernel(g::RealDataGrid, dim::Integer)
    r = fft(g.data, dim)
    for (n,s) in enumerate(eachslice(r, dims = dim))
        s .*= 2π * im * (n - 1) / size(g, dim)
    end
    return bfft!(r, dim)
end

function pdev_kernel(g::RealDataGrid{D,<:Real}, dim::Integer) where D
    r = rfft(g.data, dim)
    for (n,s) in enumerate(eachslice(r, dims = dim))
        s .*= 2π * im * (n - 1) / size(g, dim)
    end
    return brfft(r, size(g, dim), dim)
end

"""
    partial_derivative(g::RealDataGrid, dim::Integer)

Calculates the partial derivative of `g` along coordinates with respect to the basis vector of
dimension `dim`.
"""
function partial_derivative(g::RealDataGrid, dim::Integer)
    return RealDataGrid(pdev_kernel(g, dim), basis(g), shift(g))
end

"""
    cell_gradient(g::RealDataGrid{D}) -> RealDataGrid{D,<:SVector{D}}

Calculates the gradient of a `RealDataGrid`, which consists of vectors containing the components of
the partial derivatives along coordinates with respect to the basis vectors of the unit cell.
The gradient may need to be transformed to obtain a more standard result: for this, use the
`gradient()` function.
"""
function cell_gradient(g::RealDataGrid{D}) where D
    # Perform a single run to get the element type for the FFT
    pdev = pdev_kernel(g, 1)
    # Preallocate the array for storing the final result
    result = Array{SVector{D,eltype(pdev)},D}(undef, size(pdev))
    # Append the result 
    for n in eachindex(result)
        result[n] = pdev[n] * SUnitVector(1)
    end
    # Repeat this again for each index
    for i in 2:D
        pdev .= pdev_kernel(g, i)
        for n in eachindex(result)
            result[n] = pdev[n] * SUnitVector(i)
        end
    end
    return RealDataGrid(result, basis(g), shift(g))
end

cell_gradient(::RealDataGrid{0}) = error("Cannot take gradients of zero-dimensional arrays")

"""
    gradient(g::RealDataGrid{D}) -> RealDataGrid{D,<:SVector{D}}

Calculates the gradient of a `RealDataGrid`, which consists of vectors containing the components of
the partial derivatives along the orthonormal real space basis. The units of the output are those
of the input multiplied by inverse bohr.
"""
function gradient(g::RealDataGrid{D}) where D
    # Perform a single run to get the element type for the FFT
    pdev = pdev_kernel(g, 1)
    # Preallocate the array for storing the final result
    result = zeros(SVector{D,eltype(pdev)}, size(pdev))
    # Append the result 
    for n in eachindex(result)
        result[n] += basis(g) \ (pdev[n] * SUnitVector{D}(1))
    end
    # Repeat this again for each index
    for i in 2:D
        pdev .= pdev_kernel(g, i)
        for n in eachindex(result)
            result[n] += basis(g) \ (pdev[n] * SUnitVector{D}(i))
        end
    end
    return RealDataGrid(result, basis(g), shift(g))
end

gradient(::RealDataGrid{0}) = error("Cannot take gradients of zero-dimensional arrays")

"""
    directional_derivative(g::RealDataGrid{D}, v::AbstractVector) -> RealDataGrid{D}

Calculates the directional derivative of a `RealDataGrid` along a vector `v` in Cartesian
coordinates. Note that `v` does not need to be normalized, nor is it automatically normalized; 
an unnormalized input will scale the results accordingly.
"""
function directional_derivative(g::RealDataGrid{D}, v::AbstractVector) where D
    # TODO: figure out if we can do this without allocating the gradient array
    return RealDataGrid([dot(v,x) for x in gradient(g)], basis(g), shift(g))
end

#=
"""
    reinterpolate(g::RealDataGrid, dims::Integer...) -> 

Performs Fourier interpolation of `g`, returning a new grid with size equal to `dims`.
"""
function reinterpolate(g::RealDataGrid{D,T}, dims::NTuple{D,<:Integer}) where T<:Number
    fftgrid = FFTW.fft(g.data)
    truncated = fftgrid[FFTBins(dims)]
end

reinterpolate(g::RealDataGrid{D}, dims::Vararg{<:Integer,D}) = reinterpolate(g, dims)
=#

function remove_shift(g::RealDataGrid)
    # If this is an integer, removing the shift is easy
    offset_float = size(g) .* shift(g)
    offset = round(Int, test_values)
    if isapprox(offset_float, offset, atol = sqrt(eps(eltype(g))))
        # Just circular shift the backing array
        return DataGrid(circshift(g.data, -offset), basis(g))
    else
        # This will need some Fourier transform magic...
        error("Shift removal for non-trivial cases has not been implemented yet")
    end
end
