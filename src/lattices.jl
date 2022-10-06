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
    det(M) < 0 && @warn "cell vectors form a left-handed coordinate system."
    return nothing
end

# TODO: can we use views on vectors of vectors to get a matrix?
function lattice_sanity_check(vs::AbstractVector{<:AbstractVector{<:Real}}) 
    return lattice_sanity_check(hcat(vs...))
end

# TODO: deprecate in favor of `RealBasis` and `ReciprocalBasis`
"""
    BasisVectors{D} <: AbstractBasis{D}

Collection of `D` basis vectors spanning `D`-dimensional space.

Internally, a `BasisVectors{D}` is represented as an `SVector{D,SVector{D,Float64}}`. This avoids a
length declaration that would be required by an `SMatrix`.

`BasisVectors` can be indexed similarly to a matrix. However, there is one major difference: 
indexing with a single value returns an `SVector{D,Float64}`, not a value as would happen with most
subtypes of `AbstractMatrix`.
"""
struct BasisVectors{D} <: AbstractBasis{D}
    vs::SVector{D,SVector{D,Float64}}
    function BasisVectors(vs::StaticVector{D,<:StaticVector{D,<:Real}}) where D
        lattice_sanity_check(vs)
        return new{D}(vs)
    end
    function BasisVectors{D}(vs::AbstractVector{<:AbstractVector{<:Real}}) where D
        lattice_sanity_check(vs)
        return new(vs)
    end
end

BasisVectors(vs::Vararg{AbstractVector{<:Real}, D}) where D = BasisVectors{D}(hcat(vs...))
BasisVectors{D}(vs::AbstractVector{<:AbstractVector}) where D = BasisVectors{D}(hcat(vs...))

# TODO: might this be better off as an inner constructor?
function BasisVectors(M::AbstractMatrix{<:Real})
    return BasisVectors(SVector{D,SVector{D,Float64}}(M[:,n] for n in 1:D))
end

#---New RealBasis and ReciprocalBasis types-------------------------------------------------------#
"""
    RealBasis{D} <: AbstractBasis{D}

A set of real space basis vectors, assumed to be in angstroms.
"""
struct RealBasis{D} <: AbstractBasis{D}
    vs::SVector{D,SVector{D,Float64}}
    function RealBasis(vs::StaticVector{D,<:StaticVector{D,<:Real}}) where D
        lattice_sanity_check(vs)
        return new{D}(vs)
    end
    function RealBasis{D}(vs::AbstractVector{<:AbstractVector{<:Real}}) where D
        lattice_sanity_check(vs)
        return new(vs)
    end
end

"""
    ReciprocalBasis{D} <: AbstractBasis{D}

A set of reciprocal space basis vectors, assumed to be in inverse angstroms.
"""
struct ReciprocalBasis{D} <: AbstractBasis{D}
    vs::SVector{D,SVector{D,Float64}}
    function ReciprocalBasis(vs::StaticVector{D,<:StaticVector{D,<:Real}}) where D
        lattice_sanity_check(vs)
        return new{D}(vs)
    end
    function ReciprocalBasis{D}(vs::AbstractVector{<:AbstractVector{<:Real}}) where D
        lattice_sanity_check(vs)
        return new(vs)
    end
end

function (::Type{T})(M::StaticMatrix{D,D,<:Real}) where {T<:AbstractBasis,D}
    # Convert the matrix to a vector of vectors
    vs = SVector{D,SVector{D,Float64}}(M[:,n] for n in 1:D)
    # Call the inner constructor
    return T(vs)
end

function (::Type{T})(M::AbstractMatrix{<:Real}) where {T<:AbstractBasis}
    # Only allow square matrices
    @assert _allsame(size(M)) "Matrix is not square."
    # Convert the matrix to a vector of vectors
    D = size(M)[1]
    vs = SVector{D,SVector{D,Float64}}(M[:,n] for n in 1:D)
    # Call the inner constructor
    return T(vs)
end

"""
    convert(::Type{<:RealBasis}, b::ReciprocalBasis) -> RealBasis
    convert(::Type{<:ReciprocalBasis}, b::RealBasis) -> ReciprocalBasis

Converts between real space and reciprocal space representations of bases. Note that this includes
a factor of 2π that is used conventionally in crystallography: conversion from `RealBasis` to
`ReciprocalBasis` multiplies by 2π, and vice versa. This ensures that the dot products between
corresponding real and reciprocal space basis vectors are always 2π.
"""
function Base.convert(::Type{<:RealBasis}, b::ReciprocalBasis)
    return RealBasis(2π * inv(transpose(matrix(b))))
end

function Base.convert(::Type{<:ReciprocalBasis}, b::RealBasis)
    return ReciprocalBasis(inv(transpose(matrix(b))) / 2π)
end

# Tools to generate 2D and 3D lattices with given angles
#-------------------------------------------------------------------------------------------------#

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

# TODO: can we leverage QR or LU decomposition to do this generally?
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

# Fundamental methods for working with BasisVectors
#-------------------------------------------------------------------------------------------------#

# This should get a vector
Base.getindex(b::AbstractBasis, ind) = b.vs[ind]
# This should treat BasisVectors like matrices
Base.getindex(b::AbstractBasis, i1, i2) = b[i2][i1]

vectors(b::AbstractBasis) = b.vs
matrix(b::AbstractBasis{D}) where D = SMatrix{D,D,Float64}(b[m,n] for m in 1:D, n in 1:D)

# This is needed for broadcasting
Base.size(::AbstractBasis{D}) where D = (D,D)
Base.length(b::AbstractBasis) = length(b.vs)
Base.iterate(b::AbstractBasis, state) = iterate(b.vs, state) 
Base.iterate(b::AbstractBasis) = iterate(b.vs)

# TODO: does this make sense?
Base.convert(::Type{SMatrix{D,D,Float64}}, b::AbstractBasis{D}) where D = matrix(b)

# Construct zero basis
function Base.zero(T::Type{<:AbstractBasis{D}}) where D
    return T(zeros(SVector{D,SVector{D,Float64}}))
end

Base.zero(::T) where {T<:AbstractBasis} = zero(T)
Base.zeros(::Type{T}) where T<:AbstractBasis = zero(T)

# Mathematical function definitions for BasisVectors
#-------------------------------------------------------------------------------------------------#

# Definitions for multiplication/division by a scalar
# TODO: there's probably a more efficient way to do this
Base.:*(s::Number, b::AbstractBasis) = BasisVectors(matrix(b) * s)
Base.:*(b::AbstractBasis, s::Number) = s * b
Base.:/(b::AbstractBasis, s::Number) = BasisVectors(matrix(b) / s)

# And multiplication/division by vectors
Base.:*(b::AbstractBasis, v::AbstractVecOrMat) = matrix(b) * v
Base.:*(v::AbstractVecOrMat, b::AbstractBasis) = v * matrix(b)
Base.:\(b::AbstractBasis, v::AbstractVecOrMat) = matrix(b) \ v

# Unit cell metrics
#-------------------------------------------------------------------------------------------------#

"""
    cell_lengths(M::AbstractMatrix) -> Vector{Float64}

Returns the lengths of the constituent vectors in a matrix representing cell vectors.
"""
cell_lengths(M::AbstractMatrix) = [norm(M[:,n]) for n = 1:size(M,2)]
# TODO: perhaps this is not necessary anymore
# But removing it might break the API
cell_lengths(b::AbstractBasis{D}) where D = SVector{D}(norm(v) for v in b)

"""
    lengths(b::AbstractBasis) -> Vector{Float64}

Returns the lengths of the constituent vectors in a matrix representing cell vectors.
"""
lengths(b::AbstractBasis{D}) where D = SVector{D}(norm(v) for v in b)
# Get the vector lengths for anything that has a defined basis
lengths(x) = lengths(basis(x))

"""
    cell_volume(M::AbstractMatrix) -> Float64

Returns the volume of a unit cell defined by a matrix. This volume does not carry the sign
(negative for cells that do not follow the right hand rule).
"""
cell_volume(M::AbstractMatrix) = abs(det(M))
# TODO: perhaps this is not necessary anymore
# But removing it might break the API
cell_volume(b::AbstractBasis) = cell_volume(matrix(b))

"""
    volume(b::AbstractBasis) -> Float64

Returns the volume of a unit cell defined by a matrix. This volume does not carry the sign
(negative for cells that do not follow the right hand rule).
"""
volume(b::AbstractBasis) = cell_volume(matrix(b))
# Get the cell volume for anything that has a defined basis
volume(x) = volume(basis(x))

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

cell_angle_cos(b::AbstractBasis) = cell_angle_cos(matrix(b))

"""
    cell_angle_rad(b) -> Vector{Float64}

Returns the angles (in radians) between each pair of basis vectors.
"""
cell_angle_rad(b) = acos.(cell_angle_cos(b))

"""
    cell_angle_deg(b) -> Vector{Float64}

Returns the angles (in degrees) between each pair of basis vectors.
"""
cell_angle_deg(b) = acosd.(cell_angle_cos(b))

# Linear algebraic manipulation of basis vector specification
#-------------------------------------------------------------------------------------------------#

LinearAlgebra.isdiag(b::AbstractBasis) = isdiag(matrix(b))
LinearAlgebra.qr(b::AbstractBasis) = qr(matrix(b))

"""
    triangularize(l::T) where T<:AbstractBasis -> T

Converts a set of basis vectors to an upper triangular form using QR decomposition.
"""
function triangularize(b::T) where T<:AbstractBasis
    R = qr(b).R
    return T(R * diagm(sign.(diag(R))))
end

"""
    triangularize(l::T, supercell::AbstractMatrix{<:Integer}) where T<:AbstractBasis -> T

Converts a set of basis vectors to an upper triangular form using QR decomposition, with an 
included conversion to a larger supercell. The resulting matrix that describes the basis vectors
will have only positive values along the diagonal.

LAMMPS expects that basis vectors are given in this format.
"""
function triangularize(b::T, sc::AbstractMatrix{<:Integer}) where T<:AbstractBasis
    # Convert the matrix to upper triangular form using QR decomposition
    # Q is the orthogonal matrix, R is the upper triangular matrix (only need R)
    R = SMatrix{length(b),length(b),Float64}(qr(matrix(b) * sc).R)
    # Ensure the diagonal elements are positive
    return T(R * diagm(sign.(diag(R))))
end

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

Base.convert(::Type{<:RealLattice}, latt::AbstractLattice) = RealLattice(latt)
Base.convert(::Type{<:ReciprocalLattice}, latt::AbstractLattice) = ReciprocalLattice(latt)

"""
    Xtal.lattice_pair_generator_3D(M::AbstractMatrix; prim=false, ctr=:P)

Generates a pair of 3D lattices that are related by common crystallographic transformations.

Returns an `NTuple{2, BasisVectors{3}` with the first matrix containing the primitive basis and the
second containing the conventional basis.
"""
function lattice_pair_generator_3D(M::AbstractMatrix; prim=false, ctr=:P)
    # Without special centering, just return the pair of basis vectors
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

#= This docstring is broken due to @doc: https://github.com/JuliaLang/julia/issues/28834
    RealLattice{3}(M::AbstractMatrix{<:Real}; prim=false, ctr=:P)

Generates a 3-dimensional real space lattice given a set of conventional or primitive basis
vectors and centering information.B y default, inputs are assumed to describe a conventional cell.
=#
function RealLattice{3}(M::AbstractMatrix{<:Real}; prim=false, ctr=:P)
    return RealLattice(lattice_pair_generator_3D(M, prim=prim, ctr=ctr)...)
end

# Metrics, but for AbstractLattice subtypes
#-------------------------------------------------------------------------------------------------#

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
   angles_rad(L::AbstractLattice{D}; prim=false) -> SVector{D,Float64}

Returns the cell vector angles in radians.
"""
function angles_rad(L::AbstractLattice{D}; prim=false) where D
    prim ? M = prim(L) : M = conv(L)
    return SVector{D,Float64}(cell_angles_rad(M))
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
