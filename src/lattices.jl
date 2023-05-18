"""
    Electrum.lattice_sanity_check(M::AbstractMatrix)

Runs checks on a matrix intended to represent basis vectors of a crystal unit cell. Returns nothing,
but warns if the cell vectors form a left-handed coordinate system, and throws an `AssertionError`
if the cell vectors are not linearly independent.
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
lattice_sanity_check(v::AbstractVector{<:AbstractVector{<:Real}}) = lattice_sanity_check(hcat(v...))

#---New RealBasis and ReciprocalBasis types--------------------------------------------------------#
"""
    RealBasis{D} <: AbstractBasis{D}

A set of real space basis vectors, assumed to be in bohr.
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

A set of reciprocal space basis vectors, assumed to be in rad*bohr⁻¹.
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

DataSpace(::Type{RealBasis{D}}) where D = ByRealSpace{D}()
DataSpace(::Type{ReciprocalBasis{D}}) where D = ByReciprocalSpace{D}()

Base.show(io::IO, b::AbstractBasis) = print(io, typeof(b), '(', matrix(b), ')')

# Convert matrix input to a vector of vectors
function (T::Type{<:AbstractBasis{D}})(M::AbstractMatrix{<:Real}) where D
    size(M) === (D,D) || throw(
        DimensionMismatch("Cannot construct a $T from a matrix with dimensions " * string(size(M)))
    )
    return T(SVector{D,SVector{D,Float64}}(M[:,n] for n in 1:D))
end

# For statically typed arrays
function (T::Union{Type{RealBasis},Type{ReciprocalBasis}})(M::StaticMatrix{D,D,<:Real}) where D
    return T(SVector{D,SVector{D,Float64}}(M[:,n] for n in 1:D))
end

"""
    basis(x) -> AbstractBasis

Returns the basis associated with some data in `x`. By default, it's assumed to be accessible via
the property `:basis`.

The return value might be a `RealBasis` or a `ReciprocalBasis`, depending on the space in which data
is represented. Use `RealBasis(g)` or `ReciprocalBasis(g)` if a specific type is needed.
"""
basis(x) = x.basis
DataSpace(T::Type) = DataSpace(fieldtype(T, :basis))

vectors(b::AbstractBasis) = b.vs
matrix(b::AbstractBasis{D}) where D = SMatrix{D,D,Float64}(b[m,n] for m in 1:D, n in 1:D)
Base.convert(::Type{T}, b::AbstractBasis) where T<:AbstractMatrix = convert(T, matrix(b))

"""
    convert(::Type{<:AbstractBasis}, b::AbstractBasis) -> T

Converts between real space and reciprocal space representations of bases. Note that this includes a
factor of 2π that is used conventionally in crystallography: conversion from `RealBasis` to
`ReciprocalBasis` multiplies by 2π, and vice versa. This ensures that the dot products between
corresponding real and reciprocal space basis vectors are always 2π.
"""
Base.convert(T::Type{<:RealBasis}, b::ReciprocalBasis) = T(2π * inv(transpose(matrix(b))))
Base.convert(T::Type{<:ReciprocalBasis}, b::RealBasis) = T(transpose(2π * inv(matrix(b))))
(T::Union{Type{RealBasis},Type{ReciprocalBasis}})(b::AbstractBasis) = convert(T, b)
 
#---Tools to generate 2D and 3D lattices with given angles-----------------------------------------#
"""
    lattice2D(a::Real, b::Real, γ::Real) -> RealBasis{2}

Constructs a set of basis vectors in an `SMatrix` that correspond to a unit cell in 2D with the
same length and angle parameters (in degrees).

By default, the b-vector is oriented along y. This selection corresponds to the default orientation
chosen by `lattice3D()`.
"""
lattice2D(a::Real, b::Real, γ::Real) = RealBasis(SMatrix{2,2}(a*sind(γ), a*cosd(γ), 0, b))

# TODO: can we leverage QR or LU decomposition to do this generally?
"""
    lattice3D(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real) -> RealBasis{3}

Constructs a set of basis vectors in an `SMatrix` that correspond to a unit cell in 3D with the same
length and angle parameters (in degrees).

By default, the b-vector is oriented along y, and the a-vector is chosen to be perpendicular to
z, leaving the c-vector to freely vary. This selection allows for the most convenient orientation of
symmetry operations.
"""
function lattice3D(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    c1 = c*(cosd(β) - cosd(γ)*cosd(α))/sind(γ)
    c2 = c*cosd(α)
    M = SMatrix{3,3,Float64}(a*sind(γ), a*cosd(γ), 0,
                                     0,        b,  0,
                                c1, c2, sqrt(c^2 - (c1^2 + c2^2)))
    return RealBasis(M)
end

#---Fundamental methods for working with basis vectors---------------------------------------------#
# This should get a vector
Base.getindex(b::AbstractBasis, ind) = b.vs[ind]
# This should treat a basis like a matrix
Base.getindex(b::AbstractBasis, i1, i2) = b[i2][i1]

# This is needed for broadcasting
Base.size(::AbstractBasis{D}) where D = (D,D)
Base.length(b::AbstractBasis) = length(b.vs)
Base.iterate(b::AbstractBasis, state) = iterate(b.vs, state) 
Base.iterate(b::AbstractBasis) = iterate(b.vs)

# TODO: does this make sense?
Base.convert(::Type{SMatrix{D,D,Float64}}, b::AbstractBasis{D}) where D = matrix(b)

# Construct zero basis
Base.zero(T::Type{<:AbstractBasis{D}}) where D = T(zeros(SVector{D,SVector{D,Float64}}))
Base.zero(::T) where {T<:AbstractBasis} = zero(T)
Base.zeros(::Type{T}) where T<:AbstractBasis = zero(T)

#---Mathematical function definitions for basis vectors--------------------------------------------#
Base.isapprox(b1::T, b2::T; kw...) where T<:AbstractBasis = isapprox(matrix(b1), matrix(b2); kw...)

# Definitions for multiplication/division by a scalar
Base.:*(s::Number, b::T) where T<:AbstractBasis = T(matrix(b) * s)
Base.:*(b::T, s::Number) where T<:AbstractBasis = s * b
Base.:/(b::T, s::Number) where T<:AbstractBasis = T(matrix(b) / s)

# And multiplication/division by vectors
Base.:*(b::AbstractBasis, v::AbstractVecOrMat) = matrix(b) * v
Base.:*(v::AbstractVecOrMat, b::AbstractBasis) = v * matrix(b)
Base.:\(b::AbstractBasis, v::AbstractVecOrMat) = matrix(b) \ v

#---Unit cell metrics------------------------------------------------------------------------------#

"""
    Electrum.cell_lengths(M::AbstractMatrix) -> Vector{Float64}

Returns the lengths of the constituent vectors in a matrix representing cell vectors.
"""
cell_lengths(M::AbstractMatrix) = [norm(M[:,n]) for n = 1:size(M,2)]

"""
    lengths(b::AbstractBasis{D}) -> SVector{D,Float64}

Returns the lengths of the basis vectors. The units correspond to the type of the basis vectors: for
`RealBasis` the units are bohr, and for `ReciprocalBasis` the units are rad*bohr⁻¹.
"""
lengths(b::AbstractBasis{D}) where D = SVector{D}(norm(v) for v in b)

"""
    lengths(x) -> Float64

Calculate the lengths of the basis vectors associated with `x`. The units correspond to the type of
the basis vectors: for `RealBasis` the units are bohr, and for `ReciprocalBasis` the units are
rad*bohr⁻¹.
"""
lengths(x) = lengths(basis(x))

"""
    Electrum.cell_volume(M::AbstractMatrix) -> Float64

Returns the volume of a unit cell defined by a matrix. This volume does not carry the sign (negative
for cells that do not follow the right hand rule).
"""
cell_volume(M::AbstractMatrix) = abs(det(M))

"""
    volume(b::AbstractBasis) -> Float64

Returns the volume of a unit cell defined by a matrix. This volume does not carry the sign (negative
for cells that do not follow the right hand rule). The units correspond to the type of the basis 
vectors: for `RealBasis` the units are bohr³, and for `ReciprocalBasis` the units are rad³*bohr⁻³.
"""
volume(b::AbstractBasis) = cell_volume(matrix(b))

"""
    volume(x) -> Float64

Calculate the volume of the basis associated with `x`. The units correspond to the type of the basis 
vectors: for `RealBasis` the units are bohr³, and for `ReciprocalBasis` the units are rad³*bohr⁻³.
"""
volume(x) = volume(basis(x))

"""
    Electrum.generate_pairs(D::Integer) -> Vector{NTuple{2,Int}}

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
    Electrum.generate_pairs(::Type{Val{D}}) -> SVector{D*(D-1)/2, NTuple{2,Int}}

Generate pairs of integers up to `D` in ascending order in an `SVector`.
"""
function generate_pairs(::Type{Val{D}}) where D
    N = Int(D*(D-1)/2)
    return SVector{N,NTuple{2,Int}}(generate_pairs(D))
end

"""
    Electrum.cell_angles_cos(M::AbstractMatrix) -> Vector{Float64}

Generates the cosines of the unit cell angles.

The angles are generated in the correct order [α, β, γ] for 3-dimensional cells. This is achieved by
reversing the output of `Electrum.generate_pairs()`. For crystals with more spatial dimensions,
this may lead to unexpected results.
"""
function cell_angles_cos(M::AbstractMatrix)
    dimpairs = reverse(generate_pairs(size(M,1)))
    return [dot(M[:,a], M[:,b])/(norm(M[:,a])*norm(M[:,b])) for (a,b) in dimpairs]
end

"""
    angles_cos(b::AbstractBasis) -> Vector{Float64}

Generates the cosines of the unit cell angles.

The angles are generated in the correct order [α, β, γ] for 3-dimensional cells. This is achieved by
reversing the output of `Electrum.generate_pairs()`. For crystals with more spatial dimensions,
this may lead to unexpected results.
"""
angles_cos(b::AbstractBasis) = cell_angles_cos(matrix(b))

"""
    angles_rad(b) -> Vector{Float64}

Returns the angles (in radians) between each pair of basis vectors.

The angles are generated in the correct order [α, β, γ] for 3-dimensional cells. This is achieved by
reversing the output of `Electrum.generate_pairs()`. For crystals with more spatial dimensions,
this may lead to unexpected results.
"""
angles_rad(b::AbstractBasis) = acos.(angles_cos(b))

"""
    angles_deg(b) -> Vector{Float64}

Returns the angles (in degrees) between each pair of basis vectors.

The angles are generated in the correct order [α, β, γ] for 3-dimensional cells. This is achieved by
reversing the output of `Electrum.generate_pairs()`. For crystals with more spatial dimensions,
this may lead to unexpected results.
"""
angles_deg(b::AbstractBasis) = acosd.(angles_cos(b))

"""
    gram(b::AbstractBasis{D}) -> SMatrix{D,D,Float64}

Returns the Gram matrix associated with a set of basis vectors. The entries of this matrix are the
dot products associated with all possible combinations of basis vectors.
"""
gram(b::AbstractBasis) = matrix(b)' * matrix(b)

#---Linear algebraic manipulation of basis vector specification------------------------------------#

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
    triangularize(l::T, sc::AbstractMatrix{<:Integer}) where T<:AbstractBasis -> T

Converts a set of basis vectors to an upper triangular form using QR decomposition, with an included
conversion to a larger supercell. The resulting matrix that describes the basis vectors will have
only positive values along the diagonal, and therefore, is always right-handed (regardless of the
transformation matrix used).

LAMMPS expects that basis vectors are given in this format.
"""
function triangularize(b::T, sc::AbstractMatrix{<:Integer}) where T<:AbstractBasis
    # Warnings for 
    det(sc) < 0 && @warn string(
        "The transformation matrix has a negative determinant.\n",
        "However, this function always returns a right-handed basis."
    )
    det(sc) == 0 && error("supplied transformation matrix is singular")
    # Convert the matrix to upper triangular form using QR decomposition
    # Q is the orthogonal matrix, R is the upper triangular matrix (only need R)
    R = SMatrix{length(b),length(b),Float64}(qr(matrix(b) * sc).R)
    # Ensure the diagonal elements are positive
    return T(R * diagm(sign.(diag(R))))
end

# TODO: Update and document this function a bit more.
# It works, but that really isn't enough for me.
# I think this was pulled directly from the WaveTrans source code
function maxHKLindex(M::AbstractMatrix{<:Real}, ecut::Real; c = 2)
    # I think the parts below convert a set of basis vectors into their reciprocals
    #----------------------------------------------------------------------------------------------#
    cosines = cell_angles_cos(M)
    crosses = [cross(M[:,a], M[:,b]) for (a,b) in zip([2,3,1], [3,1,2])]
    triples = [dot(M[:,n], crosses[n])/(norm(crosses[n])*norm(M[:,n])) for n in 1:3]
    sines = hcat([[sqrt(1-c^2), sqrt(1-c^2) ,t] for (c,t) in zip(cosines, triples)]...)
    #----------------------------------------------------------------------------------------------#
    nbmax = sqrt(c*ecut) ./ [norm(M[:,a])*sines[a,b] for a in 1:3, b in 1:3] .+ 1
    return floor.(Int, vec(maximum(nbmax, dims=1)))
end

"""
    Electrum.maxHKLindex(b::AbstractBasis, ecut::Real; prim=true, c = CVASP)

Determines the maximum integer values of the reciprocal lattice vectors needed to store data out to
a specific energy cutoff for a 3D lattice.

By default, the energy cutoff is assumed to be in units of eV, the reciprocal lattice vector lengths
are assumed to be in angstroms, and the value of c (2m/ħ^2) is taken from VASP's default
value (which is incorrect!). Different values of c may be used for different units.

The functionality implemented here was taken from WaveTrans:
https://www.andrew.cmu.edu/user/feenstra/wavetrans/
"""
# Assume that the basis vectors are defined in reciprocal space (??)
function maxHKLindex(b::AbstractBasis{3}, ecut::Real; c = 2)
    return maxHKLindex(matrix(ReciprocalBasis(b)), ecut, c = c)
end

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
