"""
    readPOSCAR(file) -> PeriodicAtomList{3}

Reads a VASP POSCAR file. A POSCAR contains the basis vectors of the system (potentially given with
a scaling factor), the positions of all atoms as either Cartesian or reduced coordinates, and
potentially information needed to perform an ab initio MD run (currently ignored).

By default, if the provided file path is a directory, `readPOSCAR()` will read from a file named
`POSCAR` in that directory. If no path is provided, a file named `POSCAR` in the current working
directory is read.

The similar `readCONTCAR` function does much the same as `readPOSCAR`, but defaults to the file name
`CONTCAR` instead of `POSCAR`.
"""
function readPOSCAR(io::IO)
    # Skip the comment line
    readline(io)
    # Get the scaling factor
    sc = parse(Float64, strip(readline(io)))
    # Get the basis
    latt = let lns = [readline(io) for n in 1:3]
        # Get the vectors from the next three lines
        vecs = sc * [parse.(Float64, v) for v in split.(lns)]
        # Convert to RealBasis and return
        RealBasis{3}(vecs)
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
    positions = Vector{FractionalAtomPosition{3}}(undef, sum(natomtypes))
    ctr = 1
    for (n,s) in enumerate(atomnames)
        for x in 1:natomtypes[n]
            positions[ctr] = FractionalAtomPosition{3}(s, parse.(Float64, split(readline(io))))
            ctr = ctr + 1
        end
    end
    # Velocity data is ignored
    return PeriodicAtomList(latt, positions)
end

readPOSCAR(file) = open(readPOSCAR, isdir(file) ? joinpath(file, "POSCAR") : file)
readPOSCAR() = open(readPOSCAR, "POSCAR")

"""
    readCONTCAR(file) -> PeriodicAtomList{3}

Reads a VASP CONTCAR file.  A CONTCAR file is written at the end of a VASP run and contains the
atomic coordinates after the calculation completed. This is relevant for geometry optimizations or
molecular dynamics runs.

The function is broadly similar to `readPOSCAR`, but the default file names for directory arguments
or no argument is `CONTCAR` instead of `POSCAR`. For more help, see `readPOSCAR`.
"""
readCONTCAR(io::IO) = readPOSCAR(io)
readCONTCAR(file) = open(readPOSCAR, isdir(file) ? joinpath(file, "CONTCAR") : file)
readCONTCAR() = open(readPOSCAR, "CONTCAR")

"""
    writePOSCAR([file = "POSCAR"], data; names = true, dummy = false, comment)

Writes crystal data to a VASP POSCAR output. The `data` can be a `PeriodicAtomList` or an
`AbstractCrystal`. If a directory names is given instead of a file name, the data will be written to
a file named `POSCAR` in the provided directory. If no filename is provided, it is written to
`POSCAR` in the current directory.

By default, atom names are written, but this is known to break VASP 4.6. This may be overridden by
setting `names` to `false`, but this is known to cause its own incompatibility issues: it is known
that VESTA will crash if the line containing the atomic names is missing.

Dummy atoms are not are not written by default, but they may be written by setting `dummy=true`.

The first line, normally used to describe the system, may be altered by passing a printable object
to `comment`.
"""
function writePOSCAR(
    io::IO,
    list::PeriodicAtomList;
    names::Bool = true,
    dummy::Bool = false,
    comment = "Written by Electrum.jl"
)
    # Write comment line and scaling factor
    println(io, comment, "\n1")
    # Write the basis vectors
    for v in basis(list)
        println(io, (@sprintf("   % -21.16f", x) for x in v)...)
    end
    # Figure out what types of atoms there are and how many
    atom_types = atomtypes(list; dummy)
    atom_counts = [p[2] for p in atomcounts(list; dummy)]
    # If names is true, print the atom names line
    if names
        println(io, (rpad(s, 10) for s in name.(atom_types))...)
    end
    # Print the number of different types of atoms, and the `Direct` keyword
    println(io, (rpad(x, 10) for x in atom_counts)..., "\nDirect")
    for t in atom_types
        # Loop through the filtered list with one atom type
        for a in filter(a -> NamedAtom(a) == t, list.atoms)
            println(io, (@sprintf("   % -21.16f", x) for x in displacement(a))...)
        end
    end
end

writePOSCAR(io::IO, xtal::AbstractCrystal; kwargs...) = writePOSCAR4(io, xtal.atoms; kwargs...)

function writePOSCAR(file, data; kwargs...)
    f = isdir(file) ? joinpath(file, "CONTCAR") : file
    open(io -> writePOSCAR(io, data; kwargs...), f, write=true)
end

writePOSCAR(data; kwargs...) = writePOSCAR("POSCAR", data; kwargs...)

"""
    writeCONTCAR([file = "CONTCAR"], data; names = true, dummy = false, comment)

Writes crystal data to a VASP POSCAR output, but with the default filename being `CONTCAR` instead
of `POSCAR`. For more detailed help, see `writePOSCAR`.
"""
writeCONTCAR(io::IO, data; kwargs...) = writePOSCAR(io, data; kwargs...)

function writeCONTCAR(file, data; kwargs...)
    f = isdir(file) ? joinpath(file, "CONTCAR") : file
    open(io -> writePOSCAR(io, data; kwargs...), f, write=true)
end

writeCONTCAR(data; kwargs...) = writePOSCAR("CONTCAR", data; kwargs...)

"""
    writePOSCAR4([file = "POSCAR"], data; dummy = false, comment)

Writes a POSCAR file readable by VASP 4 (tested on VASP 4.6). This is equivalent to calling
`writePOSCAR([file = "POSCAR"], data; names = false, dummy = false, comment)` - in other words, the
names are never written. This can solve segmentation faults that occur with VASP 4.6, but some
software that depends on these names will misbehave (notably VESTA).

For more details, see `readPOSCAR` (and `readCONTCAR`).
"""
writePOSCAR4(args...; kwargs...) = writePOSCAR(args...; kwargs..., names=false)

# Kendall got everything done before 6 PM (2022-02-01)
"""
    readWAVECAR(file) -> PlanewaveWavefunction{3,Complex{Float32}}

Reads a WAVECAR file output from a VASP 4.6 calcuation to a `PlanewaveWavefunction`.

Information about VASP WAVECAR files and much of the code was adapted from WaveTrans (originally
written in FORTRAN): https://www.andrew.cmu.edu/user/feenstra/wavetrans/

This function is limited to WAVECAR files which have an RTAG value of 45200 (meaning the data is
given as a `Complex{Float32}`) and have only a collinear magnetic field applied, like WaveTrans. It
should also be noted that the weights of the k-points are not present in the WAVECAR file, and are
set to 1 by default.

By default, the function is verbose, with output printed for every k-point parsed, due to the large
size of the wavefunction files. If this behavior is undesirable, the `quiet` keyword argument may be
set to `true`.
"""
function readWAVECAR(io::IO; quiet = false)
    # Function to increment HKL values in place 
    function incrementHKL!(hkl::AbstractVector{<:Integer}, bounds::AbstractVector{<:AbstractRange})
        # Loop through the vector indices, but in most cases we don't need them all
        for n in eachindex(hkl)
            # Increment the current vector component
            hkl[n] = (hkl[n] + 1 in bounds[n] ? hkl[n] + 1 : minimum(bounds[n]))
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
    # Reciprocal lattice vectors
    rlatt = convert(ReciprocalBasis, RealBasis{3}([read(io, Float64) for a = 1:3, b = 1:3]))
    # Get HKL coefficient bounds (as done in WaveTrans)
    hklbounds = SVector{3,UnitRange{Int}}(-g:g for g in maxHKLindex(rlatt, ecut))
    # Bare wavefunction to be filled
    wf = PlanewaveWavefunction{3,Complex{Float32}}(rlatt, nspin, nkpt, nband, hklbounds...)
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
            wf.kpoints[kp] = [read(io, Float64) for n in 1:3]
            # Get energies and occupancies
            for b in 1:nband
                # Ordering is reversed (see PlanewaveIndex above...)
                wf.energies[b, kp, s] = read(io, Float64)
                skip(io, 8)
                wf.occupancies[b, kp, s] = read(io, Float64)
            end
            quiet || @info string(
                "Read in data for k-point ", kp, "/", nkpt, " (", npw, " planewaves/band)\n",
                "Reciprocal space coordinates: ", @sprintf("[%f %f %f]", wf.kpoints[kp]...)
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
                        sumkg = [dot(wf.kpoints[kp] + hkl, rlatt[:,n]) for n in 1:3]
                        etot = _selfdot(sumkg)/CVASP
                        # Break when the G-vector energy is below ecut
                        # This may occur immediately if the k-vector already meets the criteria
                        etot < ecut ? break : incrementHKL!(hkl, hklbounds)
                    end
                    # Store the data at the HKL index
                    # Note: data is stored first by k-points, then by bands
                    wf[s, kp, b, hkl...] = pw
                    # Increment it for the next iteration
                    incrementHKL!(hkl, hklbounds)
                end
            end
        end
    end
    return wf
end

function readWAVECAR(file; quiet = false)
    open(io -> readWAVECAR(io; quiet), isdir(file) ? joinpath(file, "WAVECAR") : file)
end

readWAVECAR(; quiet = false) = readWAVECAR("WAVECAR"; quiet)

"""
    readDOSCAR(file) -> Tuple{DensityOfStates, Vector{ProjectedDensityOfStates}}

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

readDOSCAR(file) =  open(readDOSCAR, isdir(file) ? joinpath(file, "DOSCAR") : file)
readDOSCAR() = open(readDOSCAR, "DOSCAR")

"""
    readPROCAR(file) -> FatBands{3}

Reads an lm-decomposed PROCAR file from VASP and returns its data as a `FatBands{3}`.
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
    
    # Initializes containers. 
    # Some are unused for now, with potentially useful information is commented out.
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

readPROCAR(file) = open(readPROCAR, isdir(file) ? joinpath(file, "PROCAR") : file)
readPROCAR() = open(readPROCAR, "PROCAR")

"""
    get_fermi(file) -> NamedTuple{(:fermi, :alphabeta), NTuple{2,Float64}}

Reads an OUTCAR file and returns the Fermi Energy and alpha+beta value.
"""
function get_fermi(io::IO)
    readuntil(io, "E-fermi :")
    ln = split(readline(io))
    fermi = parse.(Float64, ln[1])
    alphabeta = parse.(Float64, strip(ln[5],':'))
    return (fermi = fermi, alphabeta = alphabeta)
end

get_fermi(filename) = open(get_fermi, filename)

get_fermi() = open(get_fermi, "OUTCAR")

"""
    readKPOINTS(file) -> KPointGrid{3}

Reads a KPOINTS file to get the k-point mesh. Currently, it only supports grid-generated meshes.
"""
function readKPOINTS(io::IO)
    ln = readlines(io)
    mesh = diagm(parse.(Int64,split(ln[4])))
    # Optional 5th line contains shift of k-point mesh
    if length(ln) == 5
        shift = parse.(Float64,split(ln[5]))
    else
        shift = [0.0,0.0,0.0]
    end
    kptgrid = KPointGrid{3}(mesh, shift)
    return kptgrid
end    

readKPOINTS(file) = open(readKPOINTS, isdir(file) ? joinpath(file, "KPOINTS") : file)
readKPOINTS() = open(readKPOINTS, "KPOINTS")
