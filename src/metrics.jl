"""
    distance_periodic(basis::AbstractMatrix, v1::AbstractVector, v2::AbstractVector)

Return the distance between two vectors in a given crystal basis. This function considers all 
possible neighbors given the cell boundaries, and returns the distance between the closest pair.
"""
function distance_periodic(
    basis::AbstractMatrix,
    v1::AbstractVector,
    v2::AbstractVector,
)
    # Vectors must be the same dimensionality
    @assert length(v1) == length(v2) "Vectors are of different lengths."
    # Truncate the vector components so they're between 0 and 1
    (v1t, v2t) = map(v -> v - floor.(v), (v1, v2))
    # Get the difference vector
    dred = (v1 - v2) - round.(v1t - v2t)
    # Transform it to Cartesian coordinates and get the norm
    return norm(basis*dred)
end

"""
    delaunay_periodic(l::AtomList) -> BitMatrix

Generates the Delaunay graph of a set of points in a periodic basis.
"""
function delaunay_periodic(l::AtomList{D}) where D
    # Initialize the adjacency matrix
    M = BitMatrix(zeros(Bool, length(points), length(points)))
    # Throw an exception if there is no basis given
    basis(l) == zeros(SMatrix{D,D,Float64}) && error("No basis vectors were provided.")

    return M
end

"""
    packing_efficiency(xtal::AbstractCrystal)

Determines the packing efficiency of a crystal.
"""
function packing_efficiency(xtal::Crystal)
    
end