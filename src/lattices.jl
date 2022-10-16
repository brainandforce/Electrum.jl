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
    return RealBasis(2π * inv(transpose(matrix(b))) )
end

function Base.convert(::Type{<:ReciprocalBasis}, b::RealBasis)
    return ReciprocalBasis(transpose(2π * inv(matrix(b))))
end

RealBasis(b::AbstractBasis) = convert(RealBasis, b)
ReciprocalBasis(b::AbstractBasis) = convert(ReciprocalBasis, b)

# Tools to generate 2D and 3D lattices with given angles
#-------------------------------------------------------------------------------------------------#

"""
    lattice2D(a::Real, b::Real, γ::Real) -> RealBasis{2}

Constructs a set of basis vectors in an `SMatrix` that correspond to a unit cell in 2D with the
same length and angle parameters (in degrees).

By default, the b-vector is oriented along y. This selection corresponds to the default orientation
chosen by `lattice3D()`.
"""
function lattice2D(a::Real, b::Real, γ::Real)
    return RealBasis(SMatrix{2,2,Float64}(a*sind(γ), a*cosd(γ), 0, b))
end

# TODO: can we leverage QR or LU decomposition to do this generally?
"""
    lattice3D(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real) -> RealBasis{3}

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
    return RealBasis(M)
end

# Fundamental methods for working with basis vectors
#-------------------------------------------------------------------------------------------------#

# This should get a vector
Base.getindex(b::AbstractBasis, ind) = b.vs[ind]
# This should treat a basis like a matrix
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

# Mathematical function definitions for basis vectors
#-------------------------------------------------------------------------------------------------#

# Approximate equality
function Base.isapprox(b1::AbstractBasis, b2::AbstractBasis; kwargs...) 
    return isapprox(matrix(b1), matrix(b2); kwargs...)
end

# Definitions for multiplication/division by a scalar
Base.:*(s::Number, b::T) where T<:AbstractBasis = T(matrix(b) * s)
Base.:*(b::T, s::Number) where T<:AbstractBasis = s * b
Base.:/(b::T, s::Number) where T<:AbstractBasis = T(matrix(b) / s)

# And multiplication/division by vectors
Base.:*(b::AbstractBasis, v::AbstractVecOrMat) = matrix(b) * v
Base.:*(v::AbstractVecOrMat, b::AbstractBasis) = v * matrix(b)
Base.:\(b::AbstractBasis, v::AbstractVecOrMat) = matrix(b) \ v

# Unit cell metrics
#-------------------------------------------------------------------------------------------------#

"""
    Xtal.cell_lengths(M::AbstractMatrix) -> Vector{Float64}

Returns the lengths of the constituent vectors in a matrix representing cell vectors.
"""
cell_lengths(M::AbstractMatrix) = [norm(M[:,n]) for n = 1:size(M,2)]

"""
    lengths(b::AbstractBasis) -> Vector{Float64}

Returns the lengths of the basis vectors.
"""
lengths(b::AbstractBasis{D}) where D = SVector{D}(norm(v) for v in b)
# Get the vector lengths for anything that has a defined basis
lengths(x) = lengths(basis(x))

"""
    Xtal.cell_volume(M::AbstractMatrix) -> Float64

Returns the volume of a unit cell defined by a matrix. This volume does not carry the sign
(negative for cells that do not follow the right hand rule).
"""
cell_volume(M::AbstractMatrix) = abs(det(M))

"""
    volume(b::AbstractBasis) -> Float64

Returns the volume of a unit cell defined by a matrix. This volume does not carry the sign
(negative for cells that do not follow the right hand rule).
"""
volume(b::AbstractBasis) = cell_volume(matrix(b))
# Get the cell volume for anything that has a defined basis
volume(x) = volume(basis(x))

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
    Xtal.cell_angles_cos(M::AbstractMatrix) -> Vector{Float64}

Generates the cosines of the unit cell angles.

The angles are generated in the correct order [α, β, γ] for 3-dimensional cells. This is achieved
by reversing the output of `Xtal.generate_pairs()`. For crystals with more spatial dimensions, this
may lead to unexpected results.
"""
function cell_angles_cos(M::AbstractMatrix)
    dimpairs = reverse(generate_pairs(size(M,1)))
    return [dot(M[:,a], M[:,b])/(norm(M[:,a])*norm(M[:,b])) for (a,b) in dimpairs]
end

"""
    angles_cos(b::AbstractBasis) -> Vector{Float64}

Generates the cosines of the unit cell angles.

The angles are generated in the correct order [α, β, γ] for 3-dimensional cells. This is achieved
by reversing the output of `Xtal.generate_pairs()`. For crystals with more spatial dimensions, this
may lead to unexpected results.
"""
angles_cos(b::AbstractBasis) = cell_angles_cos(matrix(b))

"""
    angles_rad(b) -> Vector{Float64}

Returns the angles (in radians) between each pair of basis vectors.

The angles are generated in the correct order [α, β, γ] for 3-dimensional cells. This is achieved
by reversing the output of `Xtal.generate_pairs()`. For crystals with more spatial dimensions, this
may lead to unexpected results.
"""
angles_rad(b::AbstractBasis) = acos.(angles_cos(b))

"""
    angles_deg(b) -> Vector{Float64}

Returns the angles (in degrees) between each pair of basis vectors.

The angles are generated in the correct order [α, β, γ] for 3-dimensional cells. This is achieved
by reversing the output of `Xtal.generate_pairs()`. For crystals with more spatial dimensions, this
may lead to unexpected results.
"""
angles_deg(b::AbstractBasis) = acosd.(angles_cos(b))

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

# TODO: Update and document this function a bit more.
# It works, but that really isn't enough for me.
# I think this was pulled directly from the WaveTrans source code
function maxHKLindex(M::AbstractMatrix{<:Real}, ecut::Real; c = CVASP)
    # I think the parts below convert a set of basis vectors into their reciprocals
    #--------------------------------------------------------------------------#
    cosines = cell_angles_cos(M)
    crosses = [cross(M[:,a], M[:,b]) for (a,b) in zip([2,3,1], [3,1,2])]
    triples = [dot(M[:,n], crosses[n])/(norm(crosses[n])*norm(M[:,n])) for n in 1:3]
    sines = hcat([[sqrt(1-c^2), sqrt(1-c^2) ,t] for (c,t) in zip(cosines, triples)]...)
    #--------------------------------------------------------------------------#
    nbmax = sqrt(c*ecut) ./ [norm(M[:,a])*sines[a,b] for a in 1:3, b in 1:3] .+ 1
    return floor.(Int, vec(maximum(nbmax, dims=1)))
end

"""
    Xtal.maxHKLindex(b::AbstractBasis, ecut::Real; prim=true, c = CVASP)

Determines the maximum integer values of the reciprocal lattice vectors needed to store data out
to a specific energy cutoff for a 3D lattice.

By default, the energy cutoff is assumed to be in units of eV, the reciprocal lattice vector
lengths are assumed to be in angstroms, and the value of c (2m/ħ^2) is taken from VASP's default
value (which is incorrect!). Different values of c may be used for different units.

The functionality implemented here was taken from WaveTrans:
https://www.andrew.cmu.edu/user/feenstra/wavetrans/
"""
# Assume that the basis vectors are defined in reciprocal space (??)
maxHKLindex(b::AbstractBasis{3}, ecut::Real; c = CVASP) = maxHKLindex(matrix(b), ecut, c = c)

"""
    d_spacing(b::AbstractBasis, miller::AbstractVector{<:Integer}, real=true) -> Float64

Measures the real space distance between planes in a lattice given by a Miller index. By default, 
the basis vectors are assumed to be real space basis vectors.
"""
function d_spacing(b::RealBasis, miller::AbstractVector{<:Integer})
    return 1 / norm(transpose(inv(matrix(b))) * miller)
end

function d_spacing(b::ReciprocalBasis, miller::AbstractVector{<:Integer})
    return 2π / norm(matrix(b) * miller)
end

d_spacing(x, miller::AbstractVector{<:Integer}) = d_spacing(basis(x), miller)
