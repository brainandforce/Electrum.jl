Base.has_offset_axes(g::AbstractDataGrid) = true

Base.convert(T::Type{<:AbstractArray}, g::AbstractDataGrid) = convert(T, g.grid)

"""
    basis(g::AbstractDataGrid{D}) -> AbstractBasis{D}

Returns the basis associated with an `AbstractDataGrid`.

The return value might be a `RealBasis` or a `ReciprocalBasis`, depending on the space in which data
is represented. Use `RealBasis(g)` or `ReciprocalBasis(g)` if a specific type is needed.
"""
basis(g::AbstractDataGrid) = g.basis
(T::AbstractBasis)(g::AbstractDataGrid) = convert(T, basis(g))

# Automatically define the data space if possible
# If there's no basis, it'll need to be defined manually
data_space(T::Type{<:AbstractDataGrid{D}}) = data_space(fieldtype(T, basis))

Base.size(g::AbstractDataGrid) = size(g.data)
Base.size(g::AbstractDataGrid, i) = size(g.data, i)

Base.axes(g::AbstractDataGrid) = range.(0, size(g) .- 1)
Base.axes(g::AbstractDataGrid, i::Integer) = range(0, size(g, i) - 1)

# Linear indexing should be defined automatically
Base.getindex(g::AbstractDataGrid, i...) = getindex(g.data, reinterpret_index(g, i)...)
Base.setindex!(g::AbstractDataGrid, x, i...) = setindex!(g.data, x, reinterpret_index(g, i)...)

# Indexing depends on the data space trait
function CartesianIndices(g::AbstractDataGrid)
    if data_space(g) isa ReciprocalSpaceData
        return CartesianIndices(Tuple((1:n) .- (div(n, 2) + 1) for n in size(g)))
    end
    return CartesianIndices(axes(g))
end
