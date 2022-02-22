"""
    lattice_sanity_check(vecs::AbstractMatrix) -> Nothing

Tests lattice vectors to ensure they are linearly independent and that they form a right-handed
coordinate system. Returns nothing, but throws an `AssertionError` if the basis vectors are not 
linearly independent.
"""
function lattice_sanity_check(vecs::AbstractMatrix; name::AbstractString="")
    # TODO: am I ever going to use `name`?
    isempty(name) || (name = name * " ")
    d = det(vecs)
    @assert (d != 0) name * "cell vectors are not linearly independent"
    # Warn 
    (d < 0) && @warn "cell vectors form a left-handed coordinate system."
    return nothing
end

"""
    lattice_pair_check(small::AbstractMatrix, large::AbstractMatrix; tol = 1e-8) -> Bool

Checks that the larger set of lattice vectors consists of integer linear combinations of the 
smaller set of lattice vectors. Returns nothing, but throws an `AssertionError` if the check fails.

The `tol` parameter can be used to set the threshold for how close the values from solving the 
system of equations need to be to an integer to succeed (default 1e-8).
"""
function lattice_pair_check(small::AbstractMatrix, large::AbstractMatrix; tol = 1e-8)
    # Perform the sanity checks
    lattice_sanity_check.((small, large))
    # Check that all the values are close to integers
    # TODO: does this criterion work properly for reciprocal lattices?
    @assert all(x -> x - round(Int, x) < tol, small\large) "The larger set of basis vectors is \
    not expressed in terms of integer multiples of the smaller set of basis vectors."
    return nothing
end

"""
    RealLattice{D}

Describes a real space crystal lattice with primitive and conventional basis vectors.
"""
struct RealLattice{D} <: AbstractLattice{D}
    prim::SMatrix{D,D,Float64}
    conv::SMatrix{D,D,Float64}
    # Inner constructor should take any AbstractMatrix{<:Real} as input
    #=
    function RealLattice{D}(
        prim::AbstractMatrix{<:Real},
        conv::AbstractMatrix{<:Real}
    ) where D
        new(prim,conv)
    end
    =#
end

"""
    ReciprocalLattice{D}

Describes a reciprocal space crystal lattice with primitive and conventional basis vectors.
"""
struct ReciprocalLattice{D} <: AbstractLattice{D}
    prim::SMatrix{D,D,Float64}
    conv::SMatrix{D,D,Float64}
    # Inner constructor should take any AbstractMatrix{<:Real} as input
    #=
    function ReciprocalLattice{D}(
        prim::AbstractMatrix{<:Real},
        conv::AbstractMatrix{<:Real}
    ) where D
        new(prim,conv)
    end
    =#
end

# Simple constructors for lattices
# These require primitive and conventional cell vectors to be provided

function RealLattice{D}(prim::AbstractMatrix{<:Real}, conv::AbstractMatrix{<:Real}) where D
    # Perform checks on the lattice pairs
    lattice_pair_check(prim,conv)
    return RealLattice{D}(prim,conv)
end

function ReciprocalLattice{D}(prim::AbstractMatrix{<:Real}, conv::AbstractMatrix{<:Real}) where D
    # Perform checks on the lattice pairs
    lattice_pair_check(conv,prim)
    return ReciprocalLattice{D}(prim,conv)
end

function RealLattice{D}(
    prim::AbstractVector{<:AbstractVector{<:Real}},
    conv::AbstractVector{<:AbstractVector{<:Real}},
) where D
    return RealLattice{D}(hcat(prim...), hcat(conv...))
end

# Conversions between real and reciprocal lattices
# It's critical that the `RealLattice` and `ReciprocalLattice` constructors are idempotent
# so that conversion can be completely seamless

RealLattice{D}(latt::RealLattice{D}) where D = latt
ReciprocalLattice{D}(latt::ReciprocalLattice{D}) where D = latt

"""
    ReciprocalLattice{D}(latt::RealLattice{D})

Converts a real lattice to its corresponding reciprocal lattice.
"""
function ReciprocalLattice{D}(latt::RealLattice{D}) where D
    # TODO: can we define this as an involution, even with a 2pi factor?
    invlatt(M::AbstractMatrix{Real}) = collect(transpose(2*pi*inv(M)))
    return ReciprocalLattice{D}(invlatt(latt.prim), invlatt(latt.conv))
end

"""
    RealLattice{D}(latt::ReciprocalLattice{D})

Converts a reciprocal lattice to its corresponding real lattice.
"""
function RealLattice{D}(latt::ReciprocalLattice{D}) where D
    # TODO: can we define this as an involution, even with a 2pi factor?
    invlatt(M::AbstractMatrix{Real}) = collect(transpose(inv(M)/(2*pi)))
    return RealLattice{D}(invlatt(latt.prim), invlatt(latt.conv))
end

"""
    lattice_pair_generator_3D(M::AbstractMatrix; prim=false, ctr=:P)

Generates a pair of 3D lattices that are related by common transformations used in .

Returns an `NTuple{2, SMatrix{3,3,Float64}}` with the first matrix containing the primitive basis 
and the second containing the conventional basis.
"""
function lattice_pair_generator_3D(M::AbstractMatrix; prim=false, ctr=:P)
        # Without special centering just return the pair of matrices
    ctr == :P && return (M,M)
    # Figure out which type of cell the input matrix defines
    # TODO: Try to ensure this operations favors a c axis pointing along the z direction
    if prim
        Mp = M
        Mc = inv(REDUCTION_MATRIX_3D[ctr]) * M
    else
        Mc = M
        Mp = REDUCTION_MATRIX_3D[ctr] * M
    end
    return (Mp, Mc)
end

# It appears these next two docstrings are broken!
#=
"""
    RealLattice{3}(M::AbstractMatrix{<:Real}; prim=false, ctr=:P)

Generates a 3-dimensional real space lattice given a set of conventional or primitive basis
vectors and centering information.

By default, inputs are assumed to describe a conventional cell.
"""
=#
function RealLattice{3}(M::AbstractMatrix{<:Real}; prim=false, ctr=:P)
    return RealLattice{3}(lattice_pair_generator_3D(M, prim=prim, ctr=ctr)...)
end

#=
"""
    RealLattice{3}(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real; ctr=:P)

Construct a real-space crystal lattice using the lengths of the lattice basis vectors and the 
angles between them (given in degrees).

By default, the c-axis of the conventional cell will be oriented along the z-axis of Cartesian 
space, unless the cell is judged to be monoclinic.
"""
function RealLattice{3}(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real; ctr=:P)
    # TODO: make this work at some point, probably...
end
=#