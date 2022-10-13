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

function distance_periodic(b::AbstractMatrix, a1::AtomPosition, a2::AtomPosition)
    return distance_periodic(b, coord(a1), coord(a2))
end

function distance_periodic(l::AtomList)
    # Matrix with all the distances
    M = zeros(Float64, natom(l), natom(l))
    for a in 1:natom
        # Distances on diagonal should be zero/don't need to be calculated
        for b in 1:(a-1)
            dist = distance_periodic(basis(l), coord(l[a]), coord(l[b]))
            # Assign off-diagonal elements symmetrically
            M[a,b] = dist
            M[b,a] = dist
        end
    end
    return M
end

"""
    delaunay_periodic(l::AtomList) -> BitMatrix

Generates the Delaunay graph of a set of points in a periodic basis.
"""
function delaunay_periodic(l::AtomList{D}) where D
    # Initialize the adjacency matrix
    M = BitMatrix(zeros(Bool, length(points), length(points)))
    # Throw an exception if there is no basis given
    basis(l) == zeros(RealBasis{D}) && error("No basis vectors were provided.")

    return M
end

"""
    packing_efficiency(xtal::AbstractCrystal)

Determines the packing efficiency of a crystal.
"""
function packing_efficiency(xtal::Crystal)
    
end
