"""
    DensityOfStates

Contains the total density of states information.
"""
struct DensityOfStates <: AbstractDensityOfStates
    # Fermi energy
    fermi::Float64
    # Energy at each point
    energy::Vector{Float64}
    # Number of states at each energy
    dos::Vector{Float64}
    # Integrated density of states (electron count up to that energy)
    int::Vector{Float64}
    function DensityOfStates(
        fermi::Real,
        energy::AbstractVector{<:Real},
        dos::AbstractVector{<:Real},
        int::AbstractVector{<:Real}
    )
        @assert size(energy) == size(dos) == size(int) string(
            "Number of energy and DOS entries do not match."
        )
        return new(fermi, energy, dos, int)
    end
end

# Create a DOS curve without integrated DOS data
"""
    DensityOfStates(fermi::Real, energy::AbstractVector{<:Real}, dos::AbstractVector{<:Real})

Creates a `DensityOfStates` object without integrated DOS information by integrating the state 
density.
"""
function DensityOfStates(
    fermi::Real,
    energy::AbstractVector{<:Real},
    dos::AbstractVector{<:Real}
)
    # Generate an integrated DOS vector
    # TODO: perhaps we can improve this with a trapezoidal approximation?
    int = [sum(dos[1:n]) for n in 1:length(dos)]
    return DensityOfStates(fermi, energy, dos, int)
end

# Getting an index produces a tuple of the energy, 
Base.getindex(d::DensityOfStates, ind) = (d.energy[ind], d.dos[ind], d.int[ind])

# Setting up for gaussian smearing of DOS curve
function gaussian(
    dos::DensityOfStates,
    sigma::Real,
)
    # Center gaussian within energy range
    # Define gaussian
    # Permute gaussian to prevent shifting of smear output
    ctr = dos.energy[div(length(dos.energy), 2)]
    g = exp.(-((dos.energy .- ctr)/sigma).^2)
    return circshift(g, 1 - findmax(g)[2])
end

"""
    smear(dos::DensityOfStates, sigma::Real)

Smears the DOS function by convoluting it with a normalized Gaussian function with width `sigma`.
"""
function smear(
    dos::DensityOfStates,
    sigma::Real
)
    # Convolute DOS with gaussian using Fourier transforms
    gsmear = real(ifft(fft(dos.dos) .* fft(gaussian(dos, sigma))))
    # Normalize gaussian smearing so that sum(gaussian smear) = sum(DOS)
    return DensityOfStates(fermi(dos), energies(dos), (gsmear)*(sum(dos.dos)/sum(gsmear)), dos.int)
end

"""
    ProjectedDensityOfStates

Contains projected density of states information.
"""
struct ProjectedDensityOfStates
    # Fermi energy
    fermi::Float64
    # Range of energies
    energy::Vector{Float64}
    # Matrix containing projected DOS data
    # Columns are the components, rows are the energies
    dos::Matrix{Float64}
    # Integrated DOS in the same format
#    int::Matrix{Float64}
    function ProjectedDensityOfStates(
        fermi::Real,
        energy::AbstractVector{<:Real},
        dos::AbstractMatrix{<:Real}
        #int::AbstractMatrix{<:Real}
    )
        @assert length(energy) == size(dos,2) "Number of energy and DOS entries do not match."
        return new(fermi, energy, dos) #, int)
    end
end

#="""
    ProjectedDensityOfStates(
        fermi::Real,
        energy::AbstractVector{<:Real},
        dos::AbstractVector{<:Real}
    )

Creates a `ProjectedDensityOfStates` object without integrated DOS information by integrating
the state density.
"""
function ProjectedDensityOfStates(
    fermi::Real,
    energy::AbstractVector{<:Real},
    dos::AbstractMatrix{<:Real}
)
    # Generate an integrated DOS matrix
    # TODO: perhaps we can improve this with a trapezoidal approximation?
    int = hcat(vec(sum(M[:,1:n], dims=2)) for n in 1:size(M,2))
    return ProjectedDensityOfStates(fermi, energy, dos, int)
end=#

"""
    fermi(d::AbstractDensityOfStates) -> Float64

Gets the Fermi energy from DOS data. There are no guarantees on the unit of energy used!
"""
fermi(d::AbstractDensityOfStates) = d.fermi

"""
    energies(d::AbstractDensityOfStates; usefermi=false) -> Vector{Float64}

Gets the range of energies in the dataset. If `usefermi` is set to true, the energies returned will
be adjusted such that the Fermi energy is set to zero.

There are no guarantees on the unit of energy used!
"""
energies(d::AbstractDensityOfStates; usefermi=false) = d.energy .- (usefermi * d.fermi)

"""
    nelectrons(d::DensityOfStates)

Gets the approximate number of electrons that are needed to reach the Fermi level.
"""
function nelectrons(d::DensityOfStates)
    # Get the nearest entries to the Fermi energy
    E = energies(d, usefermi=true)
    # Find the closest pair of energies to the Fermi energy
    inds = partialsortperm(abs.(E), 1:2)
    # Interpolate linearly between the nearest points
    # TODO: there's probably a nicer way to do this.
    (e1, e2) = [d[i][1] for i in inds]
    (i1, i2) = [d[i][3] for i in inds]
    m = (i2 - i1)/(e2 - e1)
    b = i1 - m*e1
    return m*fermi(d) + b
end
