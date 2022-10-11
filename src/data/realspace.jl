"""
    RealSpaceDataGrid{D,T} <: AbstractRealSpaceData{D}

A data grid defined in real space, containing data of type T.

By convention, indexing of `RealSpaceDataGrid` is zero-based. This convention is used so that
datasets where the first entry corresponds to data at the origin can be indexed with zeros.
However, `getindex()` is implemented such that the dataset may be indexed by any integer, with
modulo math used to convert to an index within the grid.
"""
struct RealSpaceDataGrid{D,T} <: AbstractRealSpaceData{D}
    # Basis vectors defining the lattice
    latt::RealBasis{D}
    # Shift of the origin from the lattice
    orig::SVector{D,Float64}
    # The actual data grid
    grid::Array{T,D}
    # Inner constructor
    function RealSpaceDataGrid(
        latt::AbstractBasis{D},
        orig::AbstractVector{<:Real},
        grid::Array{T,D}
    ) where {D,T}
        @assert length(orig) == D "Origin vector does not have the lattice dimensionality"
        # Make sure the shift values lie in (-0.5, 0.5]
        orig = orig - floor.(orig)
        return new{D,T}(latt, orig, grid)
    end
end

function RealSpaceDataGrid(latt::AbstractBasis{D}, grid::Array{T,D}) where {D,T}
    return RealSpaceDataGrid(latt, zeros(SVector{D,Float64}), grid)
end

"""
    basis(g::RealSpaceDataGrid{D,T}) -> RealBasis{D}

Gets the basis vectors of a `RealSpaceDataGrid`.
"""
basis(g::RealSpaceDataGrid) = g.latt

"""
    shift(g::RealSpaceDataGrid{D,T}) -> SVector{D,Float64}

Gets the shift of the datagrid off of the origin of the basis vectors.
"""
shift(g::RealSpaceDataGrid) = g.orig

"""
    RealSpaceDataGrid(f, g::RealSpaceDataGrid)

Applies a function `f` elementwise to the grid elements of a `RealSpaceDataGrid` and returns a new
`RealSpaceDataGrid`.
"""
function RealSpaceDataGrid(f, g::RealSpaceDataGrid)
    return RealSpaceDataGrid(basis(g), shift(g), f.(grid(g)))
end

# getindex() supports arbitrary integer indices for RealSpaceDataGrid
# By convention, it's zero based, so data at fractional coordinate [0,0,0] is indexable at [0,0,0]
function Base.getindex(g::RealSpaceDataGrid, inds...)
    # Perform modulo math to get the indices (to support wraparound)
    i = mod.(inds, size(g)) .+ 1
    return getindex(grid(g), i...)
end

function Base.setindex!(g::RealSpaceDataGrid, x, inds...)
    i = mod.(inds, size(g)) .+ 1
    setindex!(grid(g), x, i...)
end

# Iterator definitions: pass through matrix iteration
Base.iterate(g::RealSpaceDataGrid, i::Integer = 1) = iterate(grid(g), i)
# Definitions for linear and Cartesian indices
Base.LinearIndices(g::RealSpaceDataGrid) = LinearIndices(grid(g)) .- 1
Base.CartesianIndices(g::RealSpaceDataGrid) = CartesianIndices(Tuple(0:n-1 for n in size(g)))

# Fast linear indexing
Base.IndexStyle(::RealSpaceDataGrid) = IndexLinear()
Base.IndexStyle(::Type{<:RealSpaceDataGrid}) = IndexLinear()

Base.keys(g::RealSpaceDataGrid) = CartesianIndices(g)

Base.eachindex(s::IndexStyle, g::RealSpaceDataGrid) = eachindex(s, g.grid)
Base.eachindex(g::RealSpaceDataGrid) = (@inline(); eachindex(IndexStyle(g), g))

"""
    grid(g::RealSpaceDataGrid{D,T}) -> Array{T,D}

Gets the array that backs a `RealSpaceDataGrid{D,T}`, which is an `Array{T,D}`.
"""
grid(g::RealSpaceDataGrid) = g.grid

Base.size(g::RealSpaceDataGrid) = size(grid(g))
Base.length(g::RealSpaceDataGrid) = length(grid(g))

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

"""
    Xtal.grid_check(g1::RealSpaceDataGrid, g2::RealSpaceDataGrid)

Performs a check on two `RealSpaceDataGrid`s to ensure that the basis, origin shift, and grid
dimensions are the same before performing mathematical operations.
"""
function grid_check(g1::RealSpaceDataGrid, g2::RealSpaceDataGrid)
    @assert basis(g1) === basis(g2) "Grid basis vectors for each grid are not identical."
    @assert shift(g1) === shift(g2) "Grid shifts from origin are not identical."
    @assert size(grid(g1)) === size(grid(g2)) "Grid sizes are different."
    return nothing
end

function Base.:+(g1::RealSpaceDataGrid{D,T1}, g2::RealSpaceDataGrid{D,T2}) where {D,T1,T2}
    # Check that the grids are identical
    grid_check(g1, g2)
    # Add the two datagrids elementwise
    newgrid = grid(g1) + grid(g2)
    return RealSpaceDataGrid(basis(g1), shift(g1), newgrid)
end

function Base.:*(g1::RealSpaceDataGrid{D,T1}, g2::RealSpaceDataGrid{D,T2}) where {D,T1,T2}
    # Check that the grids are identical
    grid_check(g1, g2)
    # Add the two datagrids elementwise
    newgrid = grid(g1) .* grid(g2)
    return RealSpaceDataGrid(basis(g1), shift(g1), newgrid)
end

function Base.:*(s::Number, g::RealSpaceDataGrid{D,T}) where {D,T}
    newgrid = s * grid(g)
    return RealSpaceDataGrid(basis(g), shift(g), newgrid)
end

Base.:*(g::RealSpaceDataGrid, s::Number) = s * g
# Multiply everything in the datagrid by -1
Base.:-(g::RealSpaceDataGrid) = RealSpaceDataGrid(basis(g), shift(g), -grid(g))
# Subtract two datagrids
Base.:-(g1::RealSpaceDataGrid, g2::RealSpaceDataGrid) = +(g1, -g2)
# Divide a datagrid by a scalar
Base.:/(g::RealSpaceDataGrid, s::Number) = RealSpaceDataGrid(basis(g), shift(g), grid(g) / s)
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
    return sum(grid(g)) * voxelsize(g)
end

"""
    integrate(f::Function, g::RealSpaceDataGrid{D,T}) -> T

Applies the function `f` pointwise to the elements of a datagrid, then integrates the grid across
all voxels.
"""
function integrate(f::Function, g::RealSpaceDataGrid)
    return sum(f.(grid(g))) * voxelsize(g)
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
