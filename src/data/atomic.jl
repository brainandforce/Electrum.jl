"""
    AtomicData{D,T}

Data associated with individual atoms in a structure.

This is a type alias for `Dict{AtomPosition{D},T}`. Keys are `AtomPosition` entries, and the values
may be of any type.
"""
const AtomicData{D,T} = Dict{AtomPosition{D},T} where {D,T}

# Note: @computed structs cannot be documented normally
# Use an @doc after the struct, like such
@computed struct SphericalHarmonic{Lmax}
    v::NTuple{(Lmax+1)^2,Float64}
    # Default constructor without parameters takes numbers directly
    function SphericalHarmonic(x::Vararg{<:Real,N}) where N
        L = sqrt(length(x)) - 1
        if isinteger(L)
            Lmax = Int(L)
        else
            throw(ArgumentError(string("Cannot determine Lmax from number of arguments.")))
        end
        return new{Lmax}(x)
    end
    # Default constructor with parameters takes any iterator
    function SphericalHarmonic{Lmax}(v) where Lmax
        @assert length(v) == (Lmax+1)^2 "For Lmax == $Lmax, iterator have length $((Lmax+1)^2)"
        return new{Lmax}(Tuple(v))
    end
end

@doc """
    SphericalHarmonic{Lmax}

Real spherical harmonic components up to `Lmax`. This can be used to describe atomic orbitals or
projections of data onto atomic sites.
""" SphericalHarmonic

SphericalHarmonic(v::SVector{N,<:Real}) where N = SphericalHarmonic(v...)
SphericalHarmonic(t::NTuple{N,<:Real}) where N = SphericalHarmonic(t...)

# Multiplication with scalars
Base.:*(x::Real, sh::SphericalHarmonic) = SphericalHarmonic(x .* sh.v)
Base.:*(sh::SphericalHarmonic, x::Real) = SphericalHarmonic(x .* sh.v)
Base.:/(sh::SphericalHarmonic, x::Real) = SphericalHarmonic(sh.v ./ x)

"""
    dot(sh1::SphericalHarmonics, sh2::SphericalHarmonics) -> Float64

Calculates the dot product between the components of two spherical harmonics, which can be used to
measure the degree of similarity between them. Note that this does not account for differences in
rotation between the spherical harmonics.
"""
function LinearAlgebra.dot(sh1::SphericalHarmonic, sh2::SphericalHarmonic)
    # Get the lengths of the smallest set of harmonics 
    # (assume all not included are zero)
    l = div(min(sizeof.((sh1, sh2))...), sizeof(Float64))
    # Restrict dot product calculation to the minimum length
    return dot(sh1.v[1:l], sh2.v[1:l])
end

"""
    Xtal.sc_ind(l::Integer, m::Integer) -> Int

Gets the associated linear index for a pair of (l,m) values used in `SphericalHarmonic`.
"""
sc_ind(l::Integer, m::Integer) = l^2 + l + 1 + m

# TODO: finish the inverse of the above function
# sc_ind(x) = 

function Base.getindex(s::SphericalHarmonic{Lmax}, l::Integer, m::Integer) where Lmax
    abs(m) <= l || error("|m| must be less than l")
    l <= Lmax || error("l exceeds lmax ($Lmax)")
    return s.v[sc_ind(l, m)]
end
