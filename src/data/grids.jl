Base.has_offset_axes(g::AbstractDataGrid) = true

Base.convert(T::Type{<:Array}, g::AbstractDataGrid) = convert(T, g.grid)

"""
    basis(g::AbstractDataGrid{D}) -> AbstractBasis{D}

Returns the basis associated with an `AbstractDataGrid`.

The return value might be a `RealBasis` or a `ReciprocalBasis`, depending on the space in which data
is represented. Use `RealBasis(g)` or `ReciprocalBasis(g)` if a specific type is needed.
"""
basis(g::AbstractDataGrid) = g.basis
(T::Type{<:AbstractBasis})(g::AbstractDataGrid) = convert(T, basis(g))

# Automatically define the data space if possible
# If there's no basis, it'll need to be defined manually
data_space(T::Type{<:AbstractDataGrid}) = data_space(fieldtype(T, basis))

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

#---Grid similarity checks-------------------------------------------------------------------------#
"""
    Electrum.grid_check(g::AbstractDataGrid{D}...) -> nothing

Checks that the basis vectors associated with a set of `AbstractDataGrid` objects are identical. It
also performs any checks specific to the data type by calling `Electrum.grid_specific_check(g...)`.
"""
function grid_check(g::AbstractDataGrid...)
    grid_specific_check(g...)
    any(h -> !isapprox(basis(first(g)), basis(h)), g) && error("Basis vectors do not match.")
    return nothing
end

"""
    Electrum.grid_specific_check(g::AbstractDataGrid...) -> Nothing

Performs extra checks that might be needed for a specific type of data grid. For any new types that
subtype `AbstractDataGrid` and require extra checks, this method should be defined. As an example,
`HKLData` has a check to ensure that the k-points associated with the data grids are identical.

By default, it performs no checks. It should always return `nothing`.
"""
grid_specific_check(g...) = nothing
# Needed to resolve method ambiguities
grid_specific_check() = nothing
