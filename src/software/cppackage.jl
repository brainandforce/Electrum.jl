
"""
    readCPcoeff(file, Lmax::Val{L}=Val{6}()) -> SphericalComponents{L}

Reads in the spherical harmonic projection coefficients from a CPpackage2 calculation.

By default, CPpackage2 gives the coefficients for spherical harmonics up to a maximum l value of 6.
"""
function readCPcoeff(io::IO, Lmax::Val{L}=Val{6}()) where L
    # All the data should be in a Vector{Float64} with this one line
    data = parse.(Float64, [v[2] for v in split.(readlines(io))])
    @debug "$(length(data)) lines in file"
    # Number of spherical harmonic coefficients
    ncoeff = (L + 1)^2
    natom = div(length(data), ncoeff)
    @debug "ncoeff = $ncoeff, natom = $natom"
    return [SphericalHarmonic{L}(data[(n - 1)*ncoeff .+ (1:ncoeff)]) for n in 1:natom]
end

readCPcoeff(filename) = open(readCPcoeff, filename)

"""
    readCPgeo(file) -> Vector{AtomPosition{3}}

Reads the atomic positions used for a CPpackage2 calculation.
"""
function readCPgeo(io::IO)
    lns = readlines(io)
    atomnames = [v[1] for v in split.(lns)]
    positions = [parse.(Float64, v[2:4]) for v in split.(lns)]
    return AtomPosition{3}.(atomnames, positions)
end

readCPgeo(filename) = open(readCPgeo, filename)

"""
    readCPcell(file) -> RealBasis{3}

Reads the basis vectors of the unit cell used for a CPpackage2 calculation.
"""
function readCPcell(io::IO)
    matrix = hcat([parse.(Float64, v) for v in split.(readlines(io))]...)
    return RealBasis{3}(matrix)
end

readCPcell(filename) = open(readCPcell, filename)