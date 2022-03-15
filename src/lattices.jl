"""
    lattice_sanity_check(vecs::AbstractMatrix) -> Nothing

Tests lattice vectors to ensure they are linearly independent and that they form a right-handed
coordinate system. Returns nothing, but throws an `AssertionError` if the basis vectors are not 
linearly independent.
"""
function lattice_sanity_check(vecs::AbstractMatrix; name::AbstractString="")
    d = det(vecs)
    @assert (d != 0) name * " "^isempty(name) * "cell vectors are not linearly independent"
    # Warn for left-handed coordinate system
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
    function RealLattice{D}(
        prim::AbstractMatrix{<:Real},
        conv::AbstractMatrix{<:Real}
    ) where D
        # Perform checks on the lattice pairs
        lattice_pair_check(prim,conv)
        new(prim,conv)
    end
end

"""
    ReciprocalLattice{D}

Describes a reciprocal space crystal lattice with primitive and conventional basis vectors.
"""
struct ReciprocalLattice{D} <: AbstractLattice{D}
    prim::SMatrix{D,D,Float64}
    conv::SMatrix{D,D,Float64}
    # Inner constructor should take any AbstractMatrix{<:Real} as input
    function ReciprocalLattice{D}(
        prim::AbstractMatrix{<:Real},
        conv::AbstractMatrix{<:Real}
    ) where D
        # Perform checks on the lattice pairs
        lattice_pair_check(conv,prim)
        new(prim,conv)
    end
end

# Get primitive and conventional lattices
# This is the preferred way of doing so
"""
    prim(l::AbstractLattice{D}) -> SMatrix{D,D,Float64}

Returns the primitive lattice in a `RealLattice` or a `ReciprocalLattice`.
"""
prim(l::AbstractLattice{D}) where D = l.prim

"""
    conv(l::AbstractLattice{D}) -> SMatrix{D,D,Float64}

Returns the conventional lattice in a `RealLattice` or a `ReciprocalLattice`.
"""
conv(l::AbstractLattice{D}) where D = l.conv

"""
    pick(l::AbstractLattice{D}, primitive::Bool) -> SMatrix{D,D,Float64}

Returns the primitive lattice if `primitive` is true, otherwise returns the conventional lattice.
"""
pick(l::AbstractLattice{D}, primitive::Bool) where D = primitive ? prim(l) : conv(l)

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
    invlatt(M::AbstractMatrix{<:Real}) = collect(transpose(2*pi*inv(M)))
    return ReciprocalLattice{D}(invlatt(latt.prim), invlatt(latt.conv))
end

"""
    RealLattice{D}(latt::ReciprocalLattice{D})

Converts a reciprocal lattice to its corresponding real lattice.
"""
function RealLattice{D}(latt::ReciprocalLattice{D}) where D
    # TODO: can we define this as an involution, even with a 2pi factor?
    invlatt(M::AbstractMatrix{<:Real}) = collect(transpose(inv(M)/(2*pi)))
    return RealLattice{D}(invlatt(latt.prim), invlatt(latt.conv))
end

# TODO: check for type stability
RealLattice(latt::AbstractLattice{D}) where D = RealLattice{D}(latt)
ReciprocalLattice(latt::AbstractLattice{D}) where D = ReciprocalLattice{D}(latt)

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

"""
    cell_lengths(M::AbstractMatrix)

Returns the lengths of the constituent vectors in a matrix representing cell vectors.
"""
cell_lengths(M::AbstractMatrix) = [norm(M[:,n]) for n = 1:size(M,2)]

"""
    cell_volume(M::AbstractMatrix)

Returns the volume of a unit cell defined by a matrix.
"""
cell_volume(M::AbstractMatrix) = abs(det(M))

"""
    generate_pairs(D::Integer) -> Vector{NTuple{2,Int}}

Generate pairs of integers up to `D` in ascending order.
"""
function generate_pairs(D::Integer)
    out = Vector{NTuple{2,Int}}(undef, Int(D*(D-1)/2))
    c = 0
    for a = 1:D
        for b = (a+1):D
            c += 1 
            out[c] = (a,b)
        end
    end
    return out
end

"""
    generate_pairs(::Type{Val{D}}) -> SVector{D*(D-1)/2, NTuple{2,Int}}

Generate pairs of integers up to `D` in ascending order in an `SVector`.
"""
function generate_pairs(::Type{Val{D}}) where D
    N = Int(D*(D-1)/2)
    return SVector{N,NTuple{2,Int}}(generate_pairs(D))
end

"""
    cell_angle_cos(M::AbstractMatrix)

Generates the cosines of the unit cell angles.

The angles are generated in the correct order [α, β, γ] for 3-dimensional cells. This is achieved
by reversing the output of `generate_pairs()`. For crystals with more spatial dimensions, this
may lead to unexpected results.
"""
function cell_angle_cos(M::AbstractMatrix)
    dimpairs = reverse(generate_pairs(size(M,1)))
    return [dot(M[:,a], M[:,b])/(norm(M[:,a])*norm(M[:,b])) for (a,b) in dimpairs]
end

"""
    cell_angle_rad(M::AbstractMatrix)

Returns the angles (in radians) between each pair of basis vectors.
"""
cell_angle_rad(M::AbstractMatrix) = acos.(cell_angle_cos(M))
cell_angle_deg(M::AbstractMatrix) = acosd.(cell_angle_cos(M))

"""
    lengths(L::AbstractLattice{D}; prim=false) -> SVector{D,Float64}

Returns the cell vector lengths.
"""
function lengths(L::AbstractLattice{D}; prim=false) where D
    prim ? M = prim(L) : M = conv(L)
    return SVector{D,Float64}(cell_lengths(M))
end

"""
    angles_deg(L::AbstractLattice{D}; prim=false) -> SVector{D,Float64}

Returns the cell vector angles in degrees.
"""
function angles_deg(L::AbstractLattice{D}; prim=false) where D
    prim ? M = prim(L) : M = conv(L)
    return SVector{D,Float64}(cell_angles_deg(M))
end

"""
   angles_deg(L::AbstractLattice{D}; prim=false) -> SVector{D,Float64}

Returns the cell vector angles in radians.
"""
function angles_rad(L::AbstractLattice{D}; prim=false) where D
    prim ? M = prim(L) : M = conv(L)
    return SVector{D,Float64}(cell_angles_rad(M))
end

"""
    lattice2D(a::Real, b::Real, γ::Real) -> SMatrix{2,2,Float64}

Constructs a set of basis vectors in an `SMatrix` that correspond to a unit cell in 2D with the
same length and angle parameters (in degrees).

By default, the b-vector is oriented along y. This selection corresponds to the default orientation
chosen by `lattice3D()`.
"""
function lattice2D(a::Real, b::Real, γ::Real)
    return SMatrix{2,2,Float64}(a*sind(γ), a*cosd(γ), 0, b)
end

"""
    lattice3D(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real) -> SMatrix{3,3,Float64}

Constructs a set of basis vectors in an `SMatrix` that correspond to a unit cell in 3D with the
same length and angle parameters (in degrees).

By default, the b-vector is oriented along y, and the a-vector is chosen to be perpendicular to
z, leaving the c-vector to freely vary. This selection allows for the most convenient orientation
of symmetry operations.
"""
function lattice3D(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    c1 = c*(cosd(β) - cosd(γ)*cosd(α))/sind(γ)
    c2 = c*cosd(α)
    return SMatrix{3,3,Float64}(a*sind(γ), a*cosd(γ), 0,
                                        0,        b,  0,
                                c1, c2, sqrt(c^2 - (c1^2 + c2^2)))
end

function maxHKLindex(M::AbstractMatrix{<:Real}, ecut::Real; c = CVASP)
    cosines = cell_angle_cos(M)
    crosses = [cross(M[:,a], M[:,b]) for (a,b) in zip([2,3,1], [3,1,2])]
    triples = [dot(M[:,n], crosses[n])/(norm(crosses[n])*norm(M[:,n])) for n in 1:3]
    sines = hcat([[sqrt(1-c^2), sqrt(1-c^2) ,t] for (c,t) in zip(cosines, triples)]...)
    nbmax = sqrt(c*ecut) ./ [norm(M[:,a])*sines[a,b] for a in 1:3, b in 1:3] .+ 1
    return floor.(Int, vec(maximum(nbmax, dims=1)))
end

"""
    maxHKLindex(L::AbstractLattice, ecut::Real; prim=true, c = CVASP)

Determines the maximum integer values of the reciprocal lattice vectors needed to store data out
to a specific energy cutoff for a 3D lattice.

By default, the energy cutoff is assumed to be in units of eV, and the value of c (2m/ħ^2) is
taken from VASP's default value. This value is off by a small amount. If different units are 
needed, the value of c should be adjusted.

The functionality implemented here would not have been possible without the work of the authors
of WaveTrans, R. M. Feenstra and M. Widom.
"""
function maxHKLindex(L::AbstractLattice{3}, ecut::Real, prim=true, c = CVASP)
    M = ReciprocalLattice{3}(prim ? prim(L) : conv(L))
    return maxHKLindex(M, ecut, c=c)
end