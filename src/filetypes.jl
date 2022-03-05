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
function readXSF3D(
    io::IO;
    spgrp::Integer = 0,
    origin::AbstractVector{<:Real} = [0, 0, 0],
    ctr::Symbol = :P
)
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

# Kendall got everything done before 6 PM (2022-02-01)
"""
    readWAVECAR(io::IO; ctr=:P) -> ReciprocalWavefunction{3,Float64,Float32}

Reads a WAVECAR file output from a VASP calcuation. Information about VASP WAVECAR files and much
of the code was pulled from the WaveTrans website (originally written in FORTRAN):
https://www.andrew.cmu.edu/user/feenstra/wavetrans/

This function is limited to WAVECAR files which have an RTAG value of 45200 (meaning the data is
given as a `Complex{Float64}`) and have only a collinear magnetic field applied, like WaveTrans.
"""
function readWAVECAR(io::IO; ctr=:P)
    # Function to increment HKL values in place 
    function incrementHKL!(hkl::AbstractVector{<:Integer}, bounds::AbstractVector{<:AbstractRange})
        # Loop through the vector indices, but in most cases we don't need them all
        for n in eachindex(hkl)
            # Increment the current vector component
            hkl[n] = (hkl[n] in bounds[n] ? hkl[n] + 1 : minimum(bounds[n]))
            # Only increment the next components if the current one is zero
            hkl[n] == 0 || break
        end
    end
    # Data entry counter (for the entries in the WAVECAR)
    count = 0
    # Number of bytes per record
    nrecl = Int(read(io, Float64))
    # Number of spin components
    nspin = Int(read(io, Float64))
    # Check for the proper format
    rtag = Int(read(io, Float64))
    rtag == 45200 || error("Unsupported format: format value is " * string(rtag))
    # Jump to the next record
    count +=1; seek(io, count*nrecl)
    # Number of k-points
    nkpt = Int(read(io, Float64))
    # Number of bands
    nband = Int(read(io, Float64))
    # Energy cutoff
    ecut = read(io, Float64)
    # Real and reciprocal lattice vectors
    latt = RealLattice{3}([read(io, Float64) for a = 1:3, b = 1:3], prim=true, ctr=ctr)
    rlatt = prim(ReciprocalLattice{3}(latt))
    # Get HKL coefficient bounds (as done in WaveTrans)
    hklbounds = SVector{3,UnitRange{Int}}(-g:g for g in maxHKLindex(rlatt, ecut))
    # List of k-points
    klist = Vector{SVector{3,Float64}}(undef, nkpt)
    # Store band info in a vector (of size nkpt) containing vectors (of size nband)
    # containing tuples (containing band energy and occupancy)
    bands = [Vector{NTuple{2,Float64}}(undef, nband) for kp in 1:nkpt]
    # Plane wave coefficients
    # Vector (size nkpt) of vectors (size nband) of vectors (size npw) of ComplexF32
    # so planewave coeffs per band per k-point - that is a *lot* of data
    waves = [[zeros(HKLData{3,Complex{Float32}}, hklbounds...) for b in 1:nband] for kp in 1:nkpt]
    # Loop through the spins
    for s in 1:nspin
        # Seek to the next data
        count += 1; seek(io, count*nrecl)
        # Loop through the k-points
        for kp in 1:nkpt
            # Number of plane waves for this k-point
            npw = read(io, Float64)
            # Add the position of the k-point to the list
            klist[kp] = [read(io, Float64) for n = 1:3]
            # Get the bands associated with the k-point
            bands[kp] = [(read(io, Float64), (skip(io, 8); read(io, Float64))) for b in 1:nband]
            for b in 1:nband
                # Seek to the next entry
                count +=1; seek(io, count*nrecl)
                # Reset the HKL indices
                hkl = zeros(MVector{3,Int})
                for p in 1:npw
                    # Get the planewave component
                    pw = read(io, Complex{Float32})
                    # Increment the HKL indices
                    while true
                        # Get the energy of the vector
                        sumkg = [dot(klist[kp] + hkl, rlatt[:,n]) for n in 1:3]
                        etot = _selfdot(sumkg)/CVASP
                        # Break when the G-vector energy is below ecut
                        # This may occur immediately if the k-vector already meets the criteria
                        etot < ecut ? break : incrementHKL!(hkl, hklbounds)
                    end
                    # Store the data at the HKL index
                    waves[kp][b][hkl...] = pw
                    # Increment it for the next iteration
                    incrementHKL!(hkl, hklbounds)
                end
            end
        end
    end
    # Now call the constructors
    return ReciprocalWavefunction{3,Float32}(
        rlatt,
        BandStructure{3}(KPointList{3}(klist), [BandAtKPoint(v) for v in bands]),
        waves
    )
end

function readWAVECAR(filename::AbstractString; ctr=:P)
    open(filename) do io
        readWAVECAR(io, ctr=ctr)
    end
end

readWAVECAR() = readWAVECAR("WAVECAR")

#=
"""
    readABINITWFK(io::IO; ctr=:P) -> (type unknown!)

Reads an ABINIT WFK file (containing a wavefunction stored by k-points).
"""
function readABINITWFK(io::IO)

end
=#