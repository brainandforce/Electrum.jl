"""
    XCFunctional

An enumerated type corresponding to different kinds of exchange-correlation functionals. The 
enumeration corresponds to the values used by ABINIT.
"""
@enum XCFunctional begin
    LDA = 1
end

"""
    HGHPseudopotential

Hartwigsen-Goedecker-Hutter pseudopotentials, as defined in "Relativistic separable dual-space 
Gaussian pseudopotentials from H to Rn" (Hartwigsen, Goedecker, Hutter, 1990) and used in 
abinit norm-conserving pseudopotential calculations.

HGH pseudopotentials are relativistic and separable, and consist of a local (real space) portion
and a nonlocal (reciprocal space) portion.
"""
struct HGHPseudopotential <: AbstractPotential
    # Atomic number
    zatom::Int8
    # Total charge
    zion::Int8
    # Exchange-correlation functional to be used
    pspxc::XCFunctional
    # Maximum value of l (for spherical harmonics)
    lmax::Int8
    # Some coefficients
    c::SVector{4,Float64}
    # r coefficients
    r::Vector{Float64}
    # h and k coefficients: dimensions 3*3*(lmax+1)
    # Note that values at l need to be accessed using the index [i,j,l+1]
    h::Array{Float64,3}
    k::Array{Float64,3}
end

"""
    zatom(p::AbstractPotential)

Gets the atomic number associated with a potential.
"""
zatom(p::AbstractPotential) = p.zatom

"""
    zion(psp::AbstractPseudopotential)

Gets the effective charge associated with a pseudopotential.
"""
zion(psp::AbstractPseudopotential) = psp.zion
zion(p::AbstractPotential) = zatom(p)

"""
    AlchemicalPotential{T<:AbstractPotential}

An alchemical potential, which can be used to describe a crystal site with mixed occupancy.
"""
struct AlchemicalPotential{T<:AbstractPotential} <: AbstractPotential
    pot::Vector{T}
    coeff::Vector{Float64}
    function AlchemicalPotential(pot::AbstractVector{T}, coeff::AbstractVector{<:Real}) where T
        @assert length(pot) == length(coeff) string(
            "Potential and coefficient vectors have different lengths"
        )
        # Remove duplicate values
    end
end