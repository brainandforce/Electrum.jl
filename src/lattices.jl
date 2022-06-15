"""
    Xtal.lattice_sanity_check(M::AbstractMatrix)

Runs checks on a matrix intended to represent basis vectors of a crystal unit cell. Returns
nothing, but warns if the cell vectors form a left-handed coordinate system, and throws an 
`AssertionError` if the cell vectors are not linearly independent.
"""
@inline function lattice_sanity_check(M::AbstractMatrix{<:Real})
    # Skip this check for lattices that are zero (meaning unspecified basis)
    if iszero(M)
        return nothing
    end
    @assert det(M) != 0 string(
        "cell vectors are not linearly independent.\n",
        "Matrix contents: ", M
    )
    det(M) < 0  && @warn "cell vectors form a left-handed coordinate system."
    return nothing
end

# TODO: can we use views on vectors of vectors to get a matrix?
function lattice_sanity_check(vs::AbstractVector{<:AbstractVector{<:Real}}) 
    return lattice_sanity_check(hcat(vs...))
end

"""
    BasisVectors{D}

Collection of `D` basis vectors spanning `D`-dimensional space.

Internally, a `BasisVectors{D}` is represented as an `SVector{D,SVector{D,Float64}}`. This avoids a
length declaration that would be required by an `SMatrix`.

`BasisVectors` can be indexed similarly to a matrix. However, there is one major difference: 
indexing with a single value returns an `SVector{D,Float64}`, not a value as would happen with most
subtypes of `AbstractMatrix`.
"""
struct BasisVectors{D}
    vs::SVector{D,SVector{D,Float64}}
    # TODO: there's probably a nicer way to handle this, but whatever
    # Inner constructor for all matrices
    function BasisVectors{D}(M::AbstractMatrix{<:Real}) where D
        lattice_sanity_check(M)
        return new{D}([M[:,n] for n in 1:size(M,2)])
    end
    # And for all vectors of vectors
    function BasisVectors{D}(M::AbstractVector{AbstractVector{<:Real}}) where D
        lattice_sanity_check(M)
        return new{D}([M[:,n] for n in 1:size(M,2)])
    end
    # If SVectors are used directly, this one can be called
    function BasisVectors(vs::SVector{D,<:SVector{D,<:Real}}) where D
        lattice_sanity_check(vs)
        return new{D}(vs)
    end
    # Same for SMatrices
    function BasisVectors(M::SMatrix{D,D,<:Real}) where D
        lattice_sanity_check(M)
        return new{D}(SVector{D,SVector{D,Float64}}(M[:,n] for n in 1:D))
    end
end

BasisVectors(vs::Vararg{AbstractVector{<:Real}, D}) where D = BasisVectors{D}(hcat(vs...))
BasisVectors{D}(vs::AbstractVector{<:AbstractVector}) where D = BasisVectors{D}(hcat(vs...))

# TODO: might this be better off as an inner constructor?
function BasisVectors(M::AbstractMatrix{<:Real})
    return BasisVectors(SVector{D,SVector{D,Float64}}(M[:,n] for n in 1:D))
end

# This should get a vector
Base.getindex(b::BasisVectors, ind) = b.vs[ind]
# This should treat BasisVectors like matrices
Base.getindex(b::BasisVectors, i1, i2) = b[i2][i1]

vectors(b::BasisVectors) = b.vs
matrix(b::BasisVectors{D}) where D = SMatrix{D,D,Float64}(b[m,n] for m in 1:D, n in 1:D)

# This is needed for broadcasting
Base.length(b::BasisVectors) = length(b.vs)
Base.iterate(b::BasisVectors, state) = iterate(b.vs, state) 
Base.iterate(b::BasisVectors) = iterate(b, 1)

# TODO: does this make sense?
Base.convert(::Type{SMatrix{D,D,Float64}}, b::BasisVectors{D}) where D = matrix(b)
# TODO: do we need other kinds of conversions?

Base.zero(::Type{BasisVectors{D}}) where D = BasisVectors(zeros(SVector{D,SVector{D,Float64}}))
Base.zero(::BasisVectors{D}) where D = BasisVectors(zeros(SVector{D,SVector{D,Float64}}))
Base.zeros(::Type{BasisVectors{D}}) where D = BasisVectors(zeros(SVector{D,SVector{D,Float64}}))

# Definitions for multiplication/division by a scalar
# TODO: there's probably a more efficient way to do this
Base.:*(s::Number, b::BasisVectors) = BasisVectors(matrix(b) * s)
Base.:*(b::BasisVectors, s::Number) = s * b
Base.:/(b::BasisVectors, s::Number) = BasisVectors(matrix(b) / s)

# And multiplication/division by vectors
Base.:*(b::BasisVectors, v::AbstractVecOrMat) = matrix(b) * v
Base.:*(v::AbstractVecOrMat, b::BasisVectors) = v * matrix(b)
Base.:\(b::BasisVectors, v::AbstractVecOrMat) = matrix(b) \ v

"""
    dual(b::BasisVectors)

Generates the dual lattice defined by a set of basis vectors.

"Dual" refers to the basis generated with the inverse transpose. This does not include factors of τ
that are used crystallographically (τ = 2π).
"""
dual(b::BasisVectors) = BasisVectors(inv(transpose(matrix(b))))

"""
    RealLattice{D}

Describes a real space crystal lattice with primitive and conventional basis vectors.
"""
struct RealLattice{D} <: AbstractLattice{D}
    prim::BasisVectors{D}
    conv::BasisVectors{D}
    # Inner constructor should take any AbstractMatrix{<:Real} as input
    function RealLattice(
        prim::BasisVectors{D},
        conv::BasisVectors{D}
    ) where D
        # Perform checks on the lattice pairs
        @assert all(x -> x - round(Int, x) < TOL_DEF, matrix(prim)\matrix(conv)) string(
            "The larger set of basis vectors is not expressed in terms of integer multiples of",
            "the smaller set of basis vectors."
        )
        return new{D}(prim, conv)
    end
end

"""
    ReciprocalLattice{D}

Describes a reciprocal space crystal lattice with primitive and conventional basis vectors.
"""
struct ReciprocalLattice{D} <: AbstractLattice{D}
    prim::BasisVectors{D}
    conv::BasisVectors{D}
    # Inner constructor should take any AbstractMatrix{<:Real} as input
    function ReciprocalLattice(
        prim::BasisVectors{D},
        conv::BasisVectors{D}
    ) where D
        # Perform checks on the lattice pairs
        @assert all(x -> x - round(Int, x) < TOL_DEF, matrix(conv)\matrix(prim)) string(
            "The smaller set of basis vectors is not expressed in terms of integer reciprocals of",
            "the larger set of basis vectors."
        )
        return new{D}(prim, conv)
    end
end

function (::Type{T})(
    b::BasisVectors{D},
    transform::AbstractMatrix=diagm(ones(D))
) where {D,T<:AbstractLattice}
    return T(b, BasisVectors{D}(b*transform))
end

# Get primitive and conventional lattices
# This is the preferred way of doing so
"""
    prim(l::AbstractLattice{D}) -> BasisVectors{D}

Returns the primitive lattice in a `RealLattice` or a `ReciprocalLattice`.
"""
prim(l::AbstractLattice) = l.prim

"""
    conv(l::AbstractLattice{D}) -> BasisVectors{D}

Returns the conventional lattice in a `RealLattice` or a `ReciprocalLattice`.
"""
conv(l::AbstractLattice) = l.conv

# Define the constructor for vectors of vectors
function (::Type{T})(
    prim::AbstractVector{<:AbstractVector{<:Real}},
    conv::AbstractVector{<:AbstractVector{<:Real}},
) where {D,T<:AbstractLattice}
    return RealLattice(BasisVectors{D}(prim), BasisVectors{D}(conv))
end

# Conversions between real and reciprocal lattices
# It's critical that the `RealLattice` and `ReciprocalLattice` constructors are idempotent
# so that conversion can be completely seamless
RealLattice(latt::RealLattice) = latt
ReciprocalLattice(latt::ReciprocalLattice) = latt

"""
    ReciprocalLattice(latt::RealLattice{D})

Converts a real lattice to its corresponding reciprocal lattice.
"""
function ReciprocalLattice(latt::RealLattice)
    return ReciprocalLattice(dual(prim(latt))*2π, dual(conv(latt))*2π)
end

"""
    RealLattice(latt::ReciprocalLattice{D})

Converts a reciprocal lattice to its corresponding real lattice.
"""
function RealLattice(latt::ReciprocalLattice)
    return RealLattice(dual(prim(latt))/2π, dual(conv(latt))/2π)
end

Base.convert(::Type{RealLattice}, latt::AbstractLattice) = RealLattice(latt)
Base.convert(::Type{ReciprocalLattice}, latt::AbstractLattice) = ReciprocalLattice(latt)
#=
function convert(::Type{RealLattice{D}}, latt::AbstractLattice{D}) where D
    return RealLattice(latt)
end

function convert(::Type{ReciprocalLattice{D}}, latt::AbstractLattice{D}) where D
    return ReciprocalLattice(latt)
end
=#

"""
    Xtal.lattice_pair_generator_3D(M::AbstractMatrix; prim=false, ctr=:P)

Generates a pair of 3D lattices that are related by common crystallographic transformations.

Returns an `NTuple{2, BasisVectors{3}` with the first matrix containing the primitive basis and the
second containing the conventional basis.
"""
function lattice_pair_generator_3D(M::AbstractMatrix; prim=false, ctr=:P)
    # Without special centering just return the pair of basis vectors
    ctr == :P && return BasisVectors{3}.((M,M))
    # Figure out which type of cell the input matrix defines
    # TODO: Try to ensure this operations favors a c axis pointing along the z direction
    if prim
        Mp = M
        Mc = inv(REDUCTION_MATRIX_3D[ctr]) * M
    else
        Mc = M
        Mp = REDUCTION_MATRIX_3D[ctr] * M
    end
    return BasisVectors{3}.((Mp,Mc))
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
    cell_lengths(M::AbstractMatrix) -> Vector{Float64}

Returns the lengths of the constituent vectors in a matrix representing cell vectors.
"""
cell_lengths(M::AbstractMatrix) = [norm(M[:,n]) for n = 1:size(M,2)]
# TODO: perhaps this is not necessary anymore
# But removing it might break the API
cell_lengths(b::BasisVectors) = cell_lengths(matrix(b))

"""
    lengths(b::BasisVectors) -> Vector{Float64}

Returns the lengths of the constituent vectors in a matrix representing cell vectors.
"""
lengths(b::BasisVectors) = cell_lengths(matrix(b))

"""
    cell_volume(M::AbstractMatrix) -> Float64

Returns the volume of a unit cell defined by a matrix. This volume does not carry the sign
(negative for cells that do not follow the right hand rule).
"""
cell_volume(M::AbstractMatrix) = abs(det(M))
# TODO: perhaps this is not necessary anymore
# But removing it might break the API
cell_volume(b::BasisVectors) = cell_volume(matrix(b))

"""
    volume(b::BasisVectors) -> Float64

Returns the volume of a unit cell defined by a matrix. This volume does not carry the sign
(negative for cells that do not follow the right hand rule).
"""
volume(b::BasisVectors) = cell_volume(matrix(b))

"""
    volume(l::AbstractLattice; primitive=true) -> Float64

Returns the volume of a lattice. By default, the primitive cell volume is used, but the
conventional cell volume may be calculated with `primitive=false`.
"""
volume(l::AbstractLattice; primitive::Bool=true) = volume(primitive ? prim(l) : conv(l))

"""
    Xtal.generate_pairs(D::Integer) -> Vector{NTuple{2,Int}}

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
    Xtal.generate_pairs(::Type{Val{D}}) -> SVector{D*(D-1)/2, NTuple{2,Int}}

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

Returns the cell vector lengths. By default, returns the lengths of the conventional cell vectors.
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
    lattice2D(a::Real, b::Real, γ::Real) -> BasisVectors{2}

Constructs a set of basis vectors in an `SMatrix` that correspond to a unit cell in 2D with the
same length and angle parameters (in degrees).

By default, the b-vector is oriented along y. This selection corresponds to the default orientation
chosen by `lattice3D()`.
"""
function lattice2D(a::Real, b::Real, γ::Real)
    return BasisVectors(SMatrix{2,2,Float64}(a*sind(γ), a*cosd(γ), 0, b))
end

"""
    lattice3D(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real) -> BasisVectors{3}

Constructs a set of basis vectors in an `SMatrix` that correspond to a unit cell in 3D with the
same length and angle parameters (in degrees).

By default, the b-vector is oriented along y, and the a-vector is chosen to be perpendicular to
z, leaving the c-vector to freely vary. This selection allows for the most convenient orientation
of symmetry operations.
"""
function lattice3D(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    c1 = c*(cosd(β) - cosd(γ)*cosd(α))/sind(γ)
    c2 = c*cosd(α)
    M = SMatrix{3,3,Float64}(a*sind(γ), a*cosd(γ), 0,
                                     0,        b,  0,
                                c1, c2, sqrt(c^2 - (c1^2 + c2^2)))
    return BasisVectors(M)
end

# TODO: Update and document this function a bit more.
# It works, but that really isn't enough for me.
# I think this was pulled directly from the WaveTrans source code
function maxHKLindex(M::AbstractMatrix{<:Real}, ecut::Real; c = CVASP)
    # I think the parts below convert a set of basis vectors into their reciprocals
    # But we already have ways to do that with `dual(::BasisVectors)`
    #--------------------------------------------------------------------------#
    cosines = cell_angle_cos(M)
    crosses = [cross(M[:,a], M[:,b]) for (a,b) in zip([2,3,1], [3,1,2])]
    triples = [dot(M[:,n], crosses[n])/(norm(crosses[n])*norm(M[:,n])) for n in 1:3]
    sines = hcat([[sqrt(1-c^2), sqrt(1-c^2) ,t] for (c,t) in zip(cosines, triples)]...)
    #--------------------------------------------------------------------------#
    nbmax = sqrt(c*ecut) ./ [norm(M[:,a])*sines[a,b] for a in 1:3, b in 1:3] .+ 1
    return floor.(Int, vec(maximum(nbmax, dims=1)))
end

# Assume that the basis vectors are defined in reciprocal space (??)
maxHKLindex(b::BasisVectors{3}, ecut::Real; c = CVASP) = maxHKLindex(matrix(b), ecut, c = c)

"""
    Xtal.maxHKLindex(L::AbstractLattice, ecut::Real; prim=true, c = CVASP)

Determines the maximum integer values of the reciprocal lattice vectors needed to store data out
to a specific energy cutoff for a 3D lattice.

By default, the energy cutoff is assumed to be in units of eV, the reciprocal lattice vector
lengths are assumed to be in angstroms, and the value of c (2m/ħ^2) is taken from VASP's default
value (which is incorrect!). Different values of c may be used for different units.

The functionality implemented here was taken from WaveTrans:
https://www.andrew.cmu.edu/user/feenstra/wavetrans/
"""
function maxHKLindex(L::AbstractLattice{3}, ecut::Real, prim=true, c = CVASP)
    M = ReciprocalLattice{3}(prim ? prim(L) : conv(L))
    return maxHKLindex(M, ecut, c=c)
end

"""
    d_spacing(b::BasisVectors, miller::AbstractVector{<:Integer}, real=true) -> Float64

Measures the real space distance between planes in a lattice given by a Miller index. By default, 
the basis vectors are assumed to be real space basis vectors.
"""
function d_spacing(b::BasisVectors, miller::AbstractVector{<:Integer}; real::Bool=true)
    return 1 / norm((real ? dual(b) : b) * miller)
end

function d_spacing(l::AbstractLattice, miller::AbstractVector{<:Integer}; primitive::Bool=false)
    r = ReciprocalLattice(l)
    return 1 / norm((primitive ? prim(r) : conv(r)) * miller)
end

function d_spacing(x::AbstractCrystal, miller::AbstractVector{<:Integer}; primitive::Bool=false)
    return d_spacing(basis(x, primitive=primitive), miller)
end
