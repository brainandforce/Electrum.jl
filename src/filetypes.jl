"""
    readXSF3D(
        io::IO;
        spgrp::Integer = 0,
        origin::AbstractVector{<:Real} = [0, 0, 0]
        ctr::Symbol = :P
    ) -> CrystalWithDatasets{3}

Reads in an XCrysDen XSF file from an input stream and returns a `CrystalWithDatasets{3}` with all
datasets that have been included within the file.

Space group and origin information are not supplied in XSF files, but they can be supplied using 
the `spgrp` and `origin` keyword arguments. Centering information can be provided using the `ctr` 
argument, but is overridden by a space group assignment.
"""
function readXSF3D(io::IO, spgrp::Integer = 0, origin::AbstractVector{<:Real})
    # Manually iterate through the file
    iter = eachline(io)
    # Counter for debugging
    ctr = 0
    # Preallocated variables
    prim = zeros(MMatrix{3,3,Float64})
    conv = zeros(MMatrix{3,3,Float64})
    atom_basis = zeros(MMatrix{3,3,Float64})
    atom_list = SVector{3,Float64}[]
    data = Dict{String,RealSpaceDataGrid{3,Float64}}()
    # There is a break statement to get of out of this loop
    while true
        # Get the line contents
        ln = let i = iterate(iter)
            # Iteration returns nothing once it's out of lines, so break from loop
            isnothing(i) ? break : i[1]
        end
        ctr += 1
        # Check for empty lines or comment lines; advance if present
        # Comments are only allowed between blocks in a valid XSF file
        # No need to check within the loops below
        (isempty(ln) || startswith(ln, "#") || all(isspace, ln)) && continue
        # The block is supposed to be given as "BEGIN_BLOCK_DATAGRID_3D"
        # But bin2xsf seems to write "BEGIN_BLOCK_DATAGRID3D" (missing underscore)
        if contains(ln, "BEGIN_BLOCK_DATAGRID") && contains(ln, "3D")
            
        elseif contains(ln, "PRIMVEC")
            
        elseif contains(ln, "CONVVEC")

        elseif contains(ln, "PRIMCOORD")

        end
    end
    xtal = Crystal{3}(RealLattice{3}(prim, conv), spgrp, origin, atomlist, atomlist)
    return CrystalWithDatasets{3,String,RealSpaceDataGrid{3}}(xtal, data)
end

"""
    readXSF3D(filename::AbstractString)

Reads an XSF file at path `filename`.
"""
function readXSF3D(filename::AbstractString)
    return open(filename) do io
        readXSF3D(io)
    end
end

# TODO: what data type do we return a WAVECAR as?
# Kendall got everything done before 6 PM (2022-02-01)
"""
    readWAVECAR(io::IO; ctr=:P) -> (type unknown!)

Reads a WAVECAR file output from a VASP calcuation. Information about VASP WAVECAR files and much
of the code was pulled from the WaveTrans website:
https://www.andrew.cmu.edu/user/feenstra/wavetrans/

This function is limited to WAVECAR files which have an RTAG value of 45200 (meaning the data is
given as a `Complex{Float64}`) and have only a collinear magnetic field applied, like WaveTrans.
"""
function readWAVECAR(io::IO; ctr=:P)
    ctr = 0
    # Number of bytes per record
    nrecl = Int(read(io, Float64))
    # Number of spin components
    nspin = Int(read(io, Float64))
    # Check for the proper format
    rtag = Int(read(io, Float64))
    rtag == 45200 || error("Unsupported format: format value is " * string(rtag))
    # Jump to the next record
    ctr +=1; seek(io, ctr*nrecl)
    # Number of k-points
    nkpt = Int(read(io, Float64))
    # Number of bands
    nband = Int(read(io, Float64))
    # Energy cutoff
    ecut = read(io, Float64)
    # Lattice vectors
    rlatt = RealLattice{3}([read(io, Float64) for a = 1:3, b = 1:3], prim=true, ctr=ctr)
    # List of k-points
    klist = Vector{SVector{3,Float64}}(undef, nkpt)
    # Store band info in a vector (of size nkpt) containing vectors (of size nband)
    # containing tuples (containing band energy and occupancy)
    bands = Vector{Vector{NTuple{2, Float64}}}(undef, nkpt)
    # Plane wave coefficients
    # Vector (size nkpt) of vectors (size nband) of vectors (size npw)
    # so planewave coeffs per band per k-point - that is a *lot* of data
    waves = Vector{Vector{Vector{Float64}}}(undef, nkpt)
    # Loop through the spins
    for s in 1:nspin
        # Seek to the next data
        ctr += 1; seek(io, ctr*nrecl)
        # Loop through the k-points
        for k in 1:nkpt
            # Number of plane waves for this k-point
            npw = read(io, Float64)
            # Add the position of the k-point to the list
            klist[k] = [read(io, Float64) for n = 1:3]
            # Get the bands associated with the k-point
            bands[k] = [(read(io, Float64), (skip(io, 8); read(io, Float64))) for b in 1:nband]
            for b in 1:nband
                # Seek to the next entry
                ctr +=1; seek(io, ctr*nrecl)
                # Get the plane wave coefficients
                # TODO: verify the data format is Complex{Float32}
                waves[k][b] = [read(io, Complex{Float32}) for w in 1:npw]
            end
        end
    end
    # Now call the constructor
    return ReciprocalWavefunction{3,Float64}(klist, bands, waves)
end

"""
    readABINITWFK(io::IO; ctr=:P) -> (type unknown!)

Reads an ABINIT WFK file (containing a wavefunction stored by k-points).
"""
function readABINITWFK(io::IO)

end