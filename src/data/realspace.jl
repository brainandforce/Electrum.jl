"""
    RealSpaceDataGrid{D,T} <: AbstractDataGrid{D,T}

A data grid defined in real space, containing data of type T.

By convention, indexing of `RealSpaceDataGrid` is zero-based. This convention is used so that
datasets where the first entry corresponds to data at the origin can be indexed with zeros. However,
`getindex()` is implemented such that the dataset may be indexed by any integer, with modulo math
used to convert to an index within the grid.
"""
struct RealSpaceDataGrid{D,T} <: AbstractDataGrid{D,T}
    # Basis vectors defining the lattice
    basis::RealBasis{D}
    # The actual data grid
    data::Array{T,D}
    # Shift of the origin from the lattice
    orig::SVector{D,Float64}
    # Inner constructor
    function RealSpaceDataGrid(
        basis::AbstractBasis{D},
        data::Array{T,D},
        orig::AbstractVector{<:Real} = zeros(SVector{D,Float64})
    ) where {D,T}
        # Make sure the shift values lie in (-0.5, 0.5]
        return new{D,T}(basis, data, orig - round.(orig))
    end
end

"""
    RealSpaceDataGrid(f::Function, g::RealSpaceDataGrid)

Applies a function `f` elementwise to the grid elements of a `RealSpaceDataGrid` and returns a new
`RealSpaceDataGrid`.
"""
RealSpaceDataGrid(f::Function, g::RealSpaceDataGrid) = RealSpaceDataGrid(g.latt, f.(g.data), g.orig)

"""
    shift(g::RealSpaceDataGrid{D,T}) -> SVector{D,Float64}

Gets the shift of the datagrid off of the origin of the basis vectors.
"""
shift(g::RealSpaceDataGrid) = g.orig

function grid_specific_check(g::RealSpaceDataGrid...)
    any(h -> !isapprox(first(g).shift, h.shift), g) && error("Grids have incommensurate shifts")
    return nothing
end

"""
    volume(g::RealSpaceDataGrid) -> Float64

Gets the crystal volume associated with a `RealSpaceDataGrid`.

By default, units are assumed to be cubic angstroms.
"""
volume(g::RealSpaceDataGrid) = volume(basis(g))

"""
    voxelsize(g::RealSpaceDataGrid) -> Float64

Gets the size of a single voxel of a `RealSpaceDataGrid`.

By default, units are assumed to be cubic angstroms.
"""
voxelsize(g::RealSpaceDataGrid) = volume(g) / length(g)

"""
    coord(g::RealSpaceDataGrid, ind...) -> SVector{D,Float64}

Returns the Cartesian coordinate associated with a grid datum at a given index.
"""
function coord(g::RealSpaceDataGrid{D,T}, ind::AbstractVector{<:Integer}) where {D,T}
    return basis(g) * (ind ./ size(g))
end

function coord(g::RealSpaceDataGrid{D,T}, ind::Vararg{<:Integer,D}) where {D,T}
    return coord(g, SVector{D,Int}(ind))
end

function coord(g::RealSpaceDataGrid{D,T}, ind::CartesianIndex{D}) where {D,T}
    return coord(g, SVector(ind.I))
end

"""
    nearest_index(g::RealSpaceDataGrid{D,T}, coord::AbstractVector{<:Real}) -> NTuple{D,Int}

Gets the nearest integer index in a data grid associated with a Cartesian coordinate.

The value returned by this function provides an index that may lie outside of the range of indices
of the backing array. However, due to the definition of `getindex()` for `RealSpaceDataGrid`, which
takes advantage of crystal periodicity, the index is guaranteed to return a value.
"""
function nearest_index(g::RealSpaceDataGrid{D,T}, coord::AbstractVector{<:Real}) where {D,T}
    return Tuple(round.(Int, (basis(g) \ coord) .* size(g)))
end

function Base.isapprox(g1::RealSpaceDataGrid, g2::RealSpaceDataGrid; kwargs...)
    grid_check(g1, g2)
    return isapprox(g1.data, g2.data, kwargs...)
end

function Base.:+(g1::RealSpaceDataGrid{D,T1}, g2::RealSpaceDataGrid{D,T2}) where {D,T1,T2}
    # Check that the grids are identical
    grid_check(g1, g2)
    # Add the two datagrids elementwise
    newgrid = g1.data .+ g2.data
    return RealSpaceDataGrid(basis(g1), newgrid, shift(g1))
end

function Base.:*(g1::RealSpaceDataGrid{D,T1}, g2::RealSpaceDataGrid{D,T2}) where {D,T1,T2}
    # Check that the grids are identical
    grid_check(g1, g2)
    # Add the two datagrids elementwise
    newgrid = g1.data .* g2.data
    return RealSpaceDataGrid(basis(g1), newgrid, shift(g1))
end

function Base.:*(s::Number, g::RealSpaceDataGrid{D,T}) where {D,T}
    newgrid = s * g.data
    return RealSpaceDataGrid(basis(g), newgrid, shift(g))
end

Base.:*(g::RealSpaceDataGrid, s::Number) = s * g
# Multiply everything in the datagrid by -1
Base.:-(g::RealSpaceDataGrid) = RealSpaceDataGrid(basis(g), -g.data, shift(g))
# Subtract two datagrids
Base.:-(g1::RealSpaceDataGrid, g2::RealSpaceDataGrid) = +(g1, -g2)
# Divide a datagrid by a scalar
Base.:/(g::RealSpaceDataGrid, s::Number) = RealSpaceDataGrid(basis(g), g.data / s, shift(g))
# Absolute value
Base.abs(g::RealSpaceDataGrid) = RealSpaceDataGrid(abs, g)
# Squared absolute value (really common!)
Base.abs2(g::RealSpaceDataGrid) = RealSpaceDataGrid(abs2, g)
# Complex angle
Base.angle(g::RealSpaceDataGrid) = RealSpaceDataGrid(g) do z
    return angle(z) + 2pi*(imag(z) < 0)
end
# Complex conjugate
Base.conj(g::RealSpaceDataGrid) = RealSpaceDataGrid(conj, g)

"""
    integrate(g::RealSpaceDataGrid{D,T}) -> T

Performs an integration across all voxels, returning a scalar value.
"""
function integrate(g::RealSpaceDataGrid)
    return sum(g.data) * voxelsize(g)
end

"""
    integrate(f::Function, g::RealSpaceDataGrid{D,T}) -> T

Applies the function `f` pointwise to the elements of a datagrid, then integrates the grid across
all voxels.
"""
function integrate(f::Function, g::RealSpaceDataGrid)
    return sum(f.(g.data)) * voxelsize(g)
end

function d_spacing(g::RealSpaceDataGrid, miller::AbstractVector{<:Integer})
    return d_spacing(basis(g), miller)
end

#=
"""
    interpolate(g::RealSpaceDataGrid{D,T}, inds...)

Evaluates a datagrid at a fractional index by linearly interpolating the values between nearest
set of points.
"""
function interpolate(g::RealSpaceDataGrid{D,T}, inds::Vararg{<:Real,D}) where {D,T}
    # Get the floor
    f = floor.(inds)
    c = ceil.(inds)
    # If they're the same, it's probably because the indices were integers
    # Just index the array
    if f == c
        return g[inds...]
    end
    # Get the set of points needed to perform the interpolation
end

"""
    interpolate(g::RealSpaceDataGrid{D,T}, r::AbstractVector{<:Real})

Evaluates a datagrid at the reduced coordinate `r` using linear interpolation.
"""
function interpolate(g::RealSpaceDataGrid, r::AbstractVector{<:Real})
    inds = reduced .* size
end
=#
