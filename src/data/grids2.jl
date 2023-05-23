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

function DataGrid{D}(data::AbstractArray{T}, basis::B, shift::S = zero_shift(b)) where {D,B,S,T}
    return DataGrid{D,B,S,T}(data, basis, shift)
end

function DataGrid(data::AbstractArray{T,D}, basis::B, shift::S = zero_shift(b)) where {D,B,S,T}
    return DataGrid{D,B,S,T}(data, basis, shift)
end

"""
    DataGrid(f, g::DataGrid{D,B,S}) -> DataGrid{D,B,S}

Constructs a new `DataGrid` by applying function or functor `f` to each element of `g`.
"""
DataGrid(f, g::DataGrid) = DataGrid(f.(g.data), g.basis, g.shift)

#---Initialize array with zeros--------------------------------------------------------------------#

function Base.zeros(
    ::Type{DataGrid{D,B,S,T}},
    basis::LatticeBasis,
    shift::AbstractVector{<:Real},
    dimensions::Vararg{<:Integer,D}
) where {D,B,S,T}
    return T(zeros(T, dimensions), basis, shift)
end

function Base.zeros(
    ::Type{DataGrid{D,B,S,T}},
    basis::LatticeBasis,
    dimensions::Vararg{<:Integer,D}
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
Base.axes(g::DataGrid) = range.(0, size(g) .- 1)

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
    voxelsize(g::RealSpaceDataGrid) -> eltype(basis(g))

Gets the size of a single real space voxel associated with a data grid. Units are assumed to be
bohr³.
"""
voxelsize(g::DataGrid) = volume(RealBasis(g)) / length(RealBasis(g))

#---Unary mathematical operations------------------------------------------------------------------#

Base.:-(g::DataGrid) = DataGrid(-g.data, g.basis, g.shift)

Base.abs(g::DataGrid) = DataGrid(abs(g.data), g.basis, g.shift)
Base.abs2(g::DataGrid) = DataGrid(abs2(g.data), g.basis, g.shift)
Base.conj(g::DataGrid) = DataGrid(conj(g.data), g.basis, g.shift)
Base.angle(g::DataGrid) = DataGrid(angle(g.data) + 2π * (imag(g.data) < 0), g.basis, g.shift)

#---Grid similarity checks-------------------------------------------------------------------------#
"""
    Electrum.grid_specific_check(g::DataGrid...) -> Nothing

Performs extra checks that might be needed for a specific type of data grid. For any new types that
subtype `AbstractDataGrid` and require extra checks, this method should be defined. As an example,
`HKLData` has a check to ensure that the k-points associated with the data grids are identical.

By default, it performs no checks. It should always return `nothing`.
"""
grid_specific_check(g...) = nothing
# Needed to resolve method ambiguities
grid_specific_check() = nothing

"""
    Electrum.grid_check(g::DataGrid{D}...) -> nothing

Checks that the basis vectors and shift associated with `DataGrid` objects are identical. It also
performs checks specific to the data type by calling `Electrum.grid_specific_check(g...)`.
"""
function grid_check(g::DataGrid...)
    grid_specific_check(g...)
    any(h -> !isapprox(basis(first(g)), basis(h)), g) && error("Basis vectors do not match.")
    return nothing
end

#---Broadcasting-----------------------------------------------------------------------------------#

function Base.similar(g::DataGrid, T::Type, sz::Tuple{Vararg{<:Integer}})
    return DataGrid(Array{T}(undef, sz), g.basis, g.shift)
end

function Base.similar(g::DataGrid, T::Type, ax::Tuple{Vararg{<:AbstractUnitRange{<:Integer}}})
    return DataGrid(Array{T}(undef, length.(ax)), g.basis, g.shift)
end

Base.similar(g::DataGrid, T::Type) = similar(g, T, size(g))
Base.similar(g::DataGrid, sz::Tuple) = similar(g, eltype(g), sz)
Base.similar(g::DataGrid) = similar(g, eltype(g), size(g))

const DataGridStyle{D,B,S} = Broadcast.ArrayStyle{DataGrid{D,B,S}}

function Base.BroadcastStyle(::Type{<:DataGrid{D,B,S}}) where {D,B,S}
    return DataGridStyle{D,B,S}()
end

# Allow broadcasting with other arrays, but don't preserve basis/shift - return Array
Base.BroadcastStyle(::DataGridStyle{D,B,S}, x::Broadcast.AbstractArrayStyle{D}) where {D,B,S} = x
Base.BroadcastStyle(::DataGridStyle{D,B,S}, x::Broadcast.DefaultArrayStyle{D}) where {D,B,S} = x

# Promote the element types of the basis and shift
function Base.BroadcastStyle(::DataGridStyle{D,B,S}, ::DataGridStyle{D,C,T}) where {D,B,S,C,T}
    return DataGridStyle{D,promote_type(B,C),promote_type(S,T)}()
end

"""
    Electrum.get_basis_and_shift(bc::Broadcast.Broadcasted)
    Electrum.get_basis_and_shift(args::NTuple{N,<:DataGrid{D,B,S} where {D,B,S}})
        -> Tuple{promote_type(B), promote_type(S)}

From a list of `DataGrid` objects, return the shared lattice basis vectors and lattice shifts. If
the set contains inequivalent values, an `Electrum.LatticeMismatch` exception will be thrown.
"""
function get_basis_and_shift(args::NTuple{N,<:DataGrid{D,B,S} where {D,B,S}}) where N
    (ub, us) = (unique(promote((getproperty(g, p) for g in args)...)) for p in (:basis, :shift))
    isone(length(ub)) || throw(LatticeMismatch("Lattice basis vectors are not equal."))
    isone(length(us)) || throw(LatticeMismatch("Lattice shifts are not equal."))
    return (only(ub), only(us))
end

get_basis_and_shift(bc::Broadcast.Broadcasted) = get_basis_and_shift(bc.args)

function Base.similar(bc::Broadcast.Broadcasted{DataGridStyle{D,B,S}}, T::Type) where {D,B,S}
    return DataGrid(Array{T}(undef, sz), get_basis_and_shift(bc)...)
end

#---More mathematical operations-------------------------------------------------------------------#

function Base.:≈(g::DataGrid, h::DataGrid; kw...)
    return all(getfield(g, f; kw...) ≈ getfield(h, f; kw...) for f in fieldnames(DataGrid))
end

Base.:*(s, g::DataGrid) = DataGrid(s * g.data, g.basis, g.shift)
Base.:*(g::DataGrid, s) = DataGrid(g.data * s, g.basis, g.shift)

function Base.:+(g::DataGrid, h::DataGrid)
    grid_check(g, h)
    return DataGrid(g.data + h.data, g.basis, g.shift)
end

function Base.:-(g::DataGrid, h::DataGrid)
    grid_check(g, h)
    return DataGrid(g.data - h.data, g.basis, g.shift)
end

#---FFTs-------------------------------------------------------------------------------------------#
"""
    fft(g::RealDataGrid) -> ReciprocalDataGrid
    fft(g::ReciprocalDataGrid) -> RealDataGrid

Performs a fast Fourier transform on the data in a `DataGrid`.
"""
FFTW.fft(g::DataGrid) = DataGrid(fft(g.data) * (det(basis(g)) / length(g)), inv(basis(g)))

#---Operations specific to `RealDataGrid`----------------------------------------------------------#
"""
    integrate([f::Function = identity], g::RealDataGrid{D,T}) -> T

Applies the function `f` elementwise to a datagrid, then integrates the grid across all voxels,
accounting for the voxel volumes.
"""
integrate(f, g::RealDataGrid) = sum(f.(g.data)) * voxelsize(g)
integrate(g::RealDataGrid) = sum(g.data) * voxelsize(g)

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
