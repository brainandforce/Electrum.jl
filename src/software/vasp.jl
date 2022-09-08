"""
    readPOSCAR(io::IO; ctr=:P) -> Crystal{3}

Reads a VASP POSCAR file.

A POSCAR contains the basis vectors of the system (potentially given with a scaling factor), the 
positions of all atoms as either Cartesian or reduced coordinates, and potentially information 
needed to perform an ab initio MD run.

Data is assumed to be given in a primitive cell, which may be converted to a conventional cell
using the `ctr` keyword.
"""
function readPOSCAR(io::IO; ctr=:P)
    # Skip the comment line
    readline(io)
    # Get the scaling factor
    sc = parse(Float64, strip(readline(io)))
    # Get the basis
    latt = let lns = [readline(io) for n in 1:3]
        # Get the vectors from the next three lines
        vecs = [parse.(Float64, v) for v in split.(lns)]
        # Return them as a lattice with desired centering
        RealLattice(BasisVectors{3}(vecs), REDUCTION_MATRIX_3D[ctr])
    end
    # Check the next line for a letter.
    # Newer versions of VASP seem to use the atom type as a line before the
    # number of each type of atom and coordintes, but version 4.6 (which we 
    # use) does not...
    ln = readline(io)
    if isletter(first(strip(ln)))
        @debug "A letter was found in the following line:\n" * ln
        atomnames = split(ln)
        ln = readline(io)
    else
        @warn string(
            "This POSCAR does not contain atomic identity information.\n",
            "Dummy atoms will be used as placeholders, with names such as X1, X2, etc."
        )
        atomnames = String[]
    end
    # Get the number of each type of atom
    natomtypes = parse.(Int, split(ln))
    if isempty(atomnames)
        atomnames = [string("X", n) for n in 1:length(natomtypes)]
    end
    # Find the "Direct" keyword
    while true
        contains("Direct", readline(io)) && break
    end
    # Generate the list of atoms
    positions = Vector{AtomPosition{3}}(undef, sum(natomtypes))
    ctr = 1
    for (n,s) in enumerate(atomnames)
        for x in 1:natomtypes[n]
            positions[ctr] = AtomPosition{3}(s, parse.(Float64, split(readline(io))))
            ctr = ctr + 1
        end
    end
    list = AtomList(latt, positions, prim=true)
    # Skip out on velocities for now
    return Crystal(latt, 1, [0, 0, 0], list, list)
end

readPOSCAR(filename::AbstractString; kwargs...) = open(readPOSCAR, filename; kwargs...)

readPOSCAR(; kwargs...) = open(readPOSCAR, "POSCAR"; kwargs...)

"""
    writePOSCAR4(io::IO, xtal::AbstractCrystal)

Writes crystal data to a VASP 4.6 POSCAR output. Dummy atoms are not written.
"""
function writePOSCAR4(io::IO, list::AtomList)
    # Write comment line and scaling factor
    println(io, "Written by Xtal.jl\n1")
    # Write the basis vectors
    for v in basis(list)
        println(io, (@sprintf("% -15.10f", x) for x in v)...)
    end
    # Figure out what types of atoms there are and how many
    atom_types = atomtypes(list)
    atom_counts = [count(a -> atomicno(a) == n, list) for n in atom_types]
    # Print the number of different types of atoms, and the `Direct` keyword
    println(io, (rpad(x, 4) for x in atom_counts)..., "\nDirect")
    for t in atom_types
        # Loop through the filtered list with one atom type
        for a in filter(a -> atomicno(a) == t, list.coord)
            println(io, (@sprintf("% -15.10f", x) for x in coord(a))...)
        end
    end
end

writePOSCAR4(io::IO, xtal::AbstractCrystal) = writePOSCAR4(io, xtal.gen)

function writePOSCAR4(filename::AbstractString, data) 
    open(filename, write=true) do io
        writePOSCAR4(io, data)
    end
end

# Kendall got everything done before 6 PM (2022-02-01)
"""
    readWAVECAR(io::IO) -> ReciprocalWavefunction{3,Float64,Float32}

Reads a WAVECAR file output from a VASP 4.6 calcuation.

Information about VASP WAVECAR files and much of the code was pulled from the WaveTrans website 
(originally written in FORTRAN): https://www.andrew.cmu.edu/user/feenstra/wavetrans/

This function is limited to WAVECAR files which have an RTAG value of 45200 (meaning the data is
given as a `Complex{Float64}`) and have only a collinear magnetic field applied, like WaveTrans. It
should also be noted that the weights of the k-points are not present in the WAVECAR file, and are
set to 1 by default.
"""
function readWAVECAR(io::IO)
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
    @debug "Record length: " * string(nrecl)
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
    latt = BasisVectors{3}([read(io, Float64) for a = 1:3, b = 1:3])
    rlatt = dual(latt) * 2pi
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
    waves = [zeros(HKLData{3,Complex{Float32}}, hklbounds...) for kp in 1:nkpt, b in 1:nband]
    # Loop through the spins
    for s in 1:nspin
        # Loop through the k-points
        for kp in 1:nkpt
            # Seek to the next data
            count += 1; seek(io, count*nrecl)
            # Number of plane waves for this k-point
            pos = position(io)
            @debug string("File pointer at ", pos, " (", count, " * ", nrecl, ")")
            npw = Int(read(io, Float64))
            # Add the position of the k-point to the list
            klist[kp] = [read(io, Float64) for n = 1:3]
            # Get the bands associated with the k-point
            bands[kp] = [(read(io, Float64), (skip(io, 8); read(io, Float64))) for b in 1:nband]
            @info string(
                "Read in data for k-point ", kp, "/", nkpt, " (", npw, " planewaves/band)\n",
                "Reciprocal space coordinates: ", @sprintf("[%f %f %f]", klist[kp]...)
            )
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
                    # Note: data is stored first by k-points, then by bands
                    waves[kp, b][hkl...] = pw
                    # Increment it for the next iteration
                    incrementHKL!(hkl, hklbounds)
                end
            end
        end
    end
    # Now call the constructors
    return ReciprocalWavefunction(rlatt, KPointList(klist), waves)
end

readWAVECAR(filename::AbstractString) = open(readWAVECAR, filename)
# Read a WAVECAR in the current directory
readWAVECAR() = readWAVECAR("WAVECAR")

"""
    readDOSCAR(io::IO) -> Tuple{DensityOfStates, Vector{ProjectedDensityOfStates}}

Reads a DOSCAR file from VASP and returns its data as a tuple containing the total and projected
density of states (if present).
"""
function readDOSCAR(io::IO)
    ln = readline(io)
    # Get number of ions and whether there is pdos data
    (nion, haspdos::Bool) = parse.(Int, split(ln))[[2,3]]
    for n in 2:6
        ln = readline(io)
    end
    # line6 = parse.(Float64, split(ln))
    # Get NEDOS and fermi energy
    nedos = parse(Int, split(ln)[3])
    fermi = parse(Float64, split(ln)[4])
    # Get the number of entries in the TDOS dataset
    ln = readline(io)
    entries = parse.(Float64, split(ln))
    tdos_raw = Matrix{Float64}(undef, length(entries), nedos)
    tdos_raw[:,1] = entries
    # Loop normally until end
    for n in 2:nedos
        ln = readline(io)
        tdos_raw[:,n] = parse.(Float64, split(ln))
    end
    # Create DensityOfStates struct for total DOS
    tdos = DensityOfStates(fermi, tdos_raw[1,:], tdos_raw[2,:], tdos_raw[3,:])
    pdos = Vector{ProjectedDensityOfStates}(undef, nion)
    if haspdos
        # Loop through all the ions
        for i in 1:nion
            readline(io)
            # Get the size of the matrix needed
            ln = readline(io)
            entries = parse.(Float64, split(ln))
            pdos_raw = Matrix{Float64}(undef, length(entries), nedos)
            pdos_raw[:,1] = entries
            # Add in the rest of the lines
            for n in 2:nedos
                ln = readline(io)
                pdos_raw[:,n] = parse.(Float64, split(ln))
            end
            # Add to the vector
            pdos[i] = ProjectedDensityOfStates(fermi, pdos_raw[1,:], pdos_raw[2:end,:])
        end
    else
        # Don't return a vector filled with a bunch of undefined objects...
        pdos = ProjectedDensityOfStates[]
    end
    return (tdos, pdos)
end

readDOSCAR(filename::AbstractString) = open(readDOSCAR, filename)

readDOSCAR() = open(readDOSCAR, "DOSCAR")

"""
    readPROCAR(io::IO) -> FatBands{3}

Reads an lm-decomposed PROCAR file from VASP and returns its data as a FatBands struct.
"""
function readPROCAR(io::IO)
    ln = readline(io)
    
    # Checks for lm-decomposed PROCAR
    if !(split(ln)[2] == "lm")
        error("Not a lm-decomposed PROCAR. Use a different PROCAR.")
    end
    has_phase = false
    if split(ln)[length(split(ln))] == "phase"
        has_phase = true
    end

    # Read header for number of kpoints, number of bands, number of ions
    (nkpt, nband, nion) = parse.(Int, split(readline(io))[[4, 8, 12]])
    
    # Initializes containers. Some are unused for now, with potentially useful information is commented out.
    # kptlist = KPointList{3}(Vector{SVector{3,Float64}}(undef, nkpt))
    # kptwt = Vector{Float64}(undef,nkpt)
    # occupancies = Matrix{Float64}(undef, nkpt, nband)
    energies = zeros(Float64, nkpt, nband)
    projband = zeros(Float64, 9, nion, nband, nkpt)
    phase = zeros(Complex{Float64}, 9, nion, nband, nkpt)
    # Begins loop
    for i in 1:nkpt
        readuntil(io, "k-point")
        ln = readline(io)
        # kptlist[i] = parse.(Float64, split(ln)[3:5])
        # kptwt[i] = parse.(Float64, split(ln)[8])
        for j in 1:nband
            readuntil(io,"band")
            ln = readline(io)
            energies[i,j] = parse(Float64, split(ln)[4])
            #occupancies[i,j] = parse(Float64, split(ln)[7])
            readuntil(io,"tot\n")
            for k in 1:nion
                ln = readline(io)
                projband[:,k,j,i] = parse.(Float64, split(ln)[2:10])
            end
            readline(io)
            if has_phase
                readline(io)
                for k in 1:nion
                    ln = readline(io)
                    re = parse.(Float64, split(ln)[2:10])
                    ln = readline(io)
                    phase[:,k,j,i] = re + parse.(Float64, split(ln)[2:10]) * im
                end
            end
        end
    end
    
    # Returns a FatBands struct
    return FatBands{3}(
        energies,
        projband,
        phase
    )
end

readPROCAR(filename::AbstractString) = open(readPROCAR, filename)

readPROCAR() = open(readPROCAR, "PROCAR")
