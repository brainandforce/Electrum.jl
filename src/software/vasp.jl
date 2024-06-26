"""
    Electrum.convert_vasp_path(path, filename)

Converts a path string (or perhaps a dedicated path type) to a filename for VASP files. Because
VASP files have predictable names, this function checks to see if the provided path is a file or
directory, then appends the expected filename to it.

The empty string is simply converted to provided filename.
"""
function convert_vasp_path(path, filename)
    isempty(path) && return(filename)
    return isdir(path) ? joinpath(path, filename) : path
end

"""
    readPOSCAR(file) -> PeriodicAtomList{3}
    readCONTCAR(file) -> PeriodicAtomList{3}

Reads a VASP POSCAR or CONTCAR file. A POSCAR contains the basis vectors of the system (potentially
given with a scaling factor), the positions of all atoms as either Cartesian or reduced coordinates,
and potentially information needed to perform an ab initio MD run (currently ignored). The similarly
formatted CONTCAR file is the geometry output generated after a VASP run.

Regardless of whether the coordinates are Cartesian or fractional, this function always returns a 
`PeriodicAtomList{3}`, and the units are converted from angstroms to bohr.

By default, if the provided file path is a directory, `readPOSCAR()` will read from a file named
`POSCAR` in that directory. If no argument is provided, a file named `POSCAR` in the current working
directory is read. `readCONTCAR()` will search for a file named `CONTCAR` in an analogous manner.
"""
function readPOSCAR(io::IO)
    # Skip the comment line
    readline(io)
    # Get the scaling factor
    sc = parse(Float64, strip(readline(io)))
    # Get the basis vectors from the next three lines and apply the scaling factor
    ln = split(join(readline(io) for _ in 1:3))
    latt = RealBasis(sc * ANG2BOHR * SMatrix{3,3}(parse.(Float64, ln)))
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
    # Create the iterator for atom names
    name_iter = Iterators.flatten(Iterators.repeated(a, b) for (a,b) in zip(atomnames, natomtypes))
    # Determine whether the coordinates are fractional (Direct) or Cartesian by looping for keyword
    # Velocity data will always be ignored
    while true
        ln = readline(io)
        if contains("direct", lowercase(ln))
            positions = Vector{FractionalAtomPosition{3}}(undef, sum(natomtypes))
            for (n,s) in enumerate(name_iter)
                ln = readline(io)
                positions[n] = FractionalAtomPosition{3}(s, parse.(Float64, split(ln)))
            end
            return PeriodicAtomList(latt, positions)
        elseif contains("cartesian", lowercase(ln))
            positions = Vector{CartesianAtomPosition{3}}(undef, sum(natomtypes))
            for (n,s) in enumerate(name_iter)
                ln = readline(io)
                # Convert the units to Bohr here
                positions[n] = CartesianAtomPosition{3}(s, ANG2BOHR * parse.(Float64, split(ln)))
            end
            return PeriodicAtomList(latt, positions)
        end
    end
    error(
        "No keyword for the atomic coordinate representation was found.\n" *
        "Valid POSCAR files should contain the keyword \"Direct\" for fractional coordinates, or " *
        "\"Cartesian\" for Cartesian coordinates."
    )
end

readPOSCAR(file) = open(readPOSCAR, convert_vasp_path(file, "POSCAR"))
readPOSCAR() = open(readPOSCAR, "POSCAR")

readCONTCAR(io::IO) = readPOSCAR(io)
readCONTCAR(file) = open(readPOSCAR, convert_vasp_path(file, "CONTCAR"))
readCONTCAR() = open(readPOSCAR, "CONTCAR")
@doc (@doc readPOSCAR) readCONTCAR

"""
    writePOSCAR([file = "POSCAR"], data; names = true, dummy = false, comment)
    writeCONTCAR([file = "CONTCAR"], data; names = true, dummy = false, comment)

Writes crystal data to a VASP POSCAR or CONTCAR output. The `data` can be a `PeriodicAtomList` or an
`AbstractCrystal`. If a directory name is given instead of a file name, the data will be written to
a file named `POSCAR` or `CONTCAR` in the provided directory, depending on the function used. If no
argument is provided, it is written to `POSCAR` or `CONTCAR` in the current directory.

By default, atom names are written, but this is known to break VASP 4.6. This may be overridden by
setting `names` to `false`, but this is known to cause its own incompatibility issues: it is known
that VESTA, for instance, will crash if the line containing the atomic names is missing.

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
    for v in eachcol(basis(list))
        println(io, (@sprintf("   % -21.16f", x * BOHR2ANG) for x in v)...)
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

writePOSCAR(io::IO, xtal::AbstractCrystal; kwargs...) = writePOSCAR(io, xtal.atoms; kwargs...)

function writePOSCAR(file, data; kwargs...)
    open(io -> writePOSCAR(io, data; kwargs...), convert_vasp_path(file, "POSCAR"), write=true)
end

writePOSCAR(data; kwargs...) = writePOSCAR("POSCAR", data; kwargs...)

writeCONTCAR(io::IO, data; kwargs...) = writePOSCAR(io, data; kwargs...)

function writeCONTCAR(file, data; kwargs...)
    open(io -> writePOSCAR(io, data; kwargs...), convert_vasp_path(file, "CONTCAR"), write=true)
end

writeCONTCAR(data; kwargs...) = writePOSCAR("CONTCAR", data; kwargs...)

@doc (@doc writePOSCAR) writeCONTCAR

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
    ct = 0
    # Number of bytes per record
    nrecl = Int(read(io, Float64))
    @debug "Record length: " * string(nrecl)
    # Number of spin components
    nspin = Int(read(io, Float64))
    # Check for the proper format
    rtag = Int(read(io, Float64))
    rtag == 45200 || error("Unsupported format: format value is " * string(rtag))
    # Jump to the next record
    ct +=1; seek(io, ct*nrecl)
    # Number of k-points
    nkpt = Int(read(io, Float64))
    # Number of bands
    nband = Int(read(io, Float64))
    # Energy cutoff
    ecut_eV = read(io, Float64)
    ecut = ecut_eV * EV2HARTREE
    # Real and reciprocal lattice vectors
    # Given in angstroms in the WAVECAR - convert to bohr for the PlanewaveWavefunction constructor
    # But use the raw data for the calculation of the G-vector bounds/array allocation
    raw_real_lattice = [read(io, Float64) for _ in 1:3, _ in 1:3]
    raw_reciprocal_lattice = 2π * inv(transpose(raw_real_lattice))
    rlatt = ReciprocalBasis{3}(raw_real_lattice / ANG2BOHR)
    # Get HKL coefficient bounds (as done in WaveTrans)
    hklbounds = SVector{3}(-g:g for g in maxHKLindex(raw_reciprocal_lattice, ecut_eV, c = CVASP))
    # Bare wavefunction to be filled
    wf = PlanewaveWavefunction{3,Complex{Float32}}(rlatt, nspin, nkpt, nband, hklbounds...)
    # Loop through the spins
    for s in 1:nspin
        # Loop through the k-points
        for kp in 1:nkpt
            # Seek to the next data
            ct += 1; seek(io, ct*nrecl)
            # Number of plane waves for this k-point
            @debug string("File pointer at ", position(io), " (", ct, " * ", nrecl, ")")
            npw = Int(read(io, Float64))
            # Add the position of the k-point to the list
            wf.kpoints[kp] = [read(io, Float64) for n in 1:3]
            # Get energies and occupancies
            for b in 1:nband
                # Ordering is reversed (see PlanewaveIndex above...)
                wf.energies[b, kp, s] = read(io, Float64) * EV2HARTREE
                skip(io, 8)
                wf.occupancies[b, kp, s] = read(io, Float64)
            end
            quiet || @info string(
                "Read in data for k-point ", kp, "/", nkpt, " (", npw, " planewaves/band)\n",
                "Reciprocal space coordinates: ", @sprintf("[%f %f %f]", wf.kpoints[kp]...)
            )
            for b in 1:nband
                # Seek to the next entry
                ct +=1; seek(io, ct*nrecl)
                # Reset the HKL indices
                hkl = zeros(MVector{3,Int})
                for p in 1:npw
                    # Get the planewave component
                    pw = read(io, Complex{Float32})
                    # Increment the HKL indices
                    while true
                        # Get the energy of the vector
                        sumkg = raw_reciprocal_lattice * (wf.kpoints[kp] + hkl)
                        etot = dot(sumkg, sumkg) / CVASP
                        # Break when the G-vector energy is below ecut
                        # This may occur immediately if the G-vector already meets the criteria
                        etot < ecut_eV ? break : incrementHKL!(hkl, hklbounds)
                    end
                    wf[s, kp, b, hkl...] = pw
                    incrementHKL!(hkl, hklbounds)
                end
            end
        end
    end
    return wf
end

function readWAVECAR(file; quiet = false)
    open(io -> readWAVECAR(io; quiet), convert_vasp_path(file, "WAVECAR"))
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
    tdos = DensityOfStates(fermi*EV2HARTREE, tdos_raw[1,:].*EV2HARTREE, tdos_raw[2,:], tdos_raw[3,:])
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
            pdos[i] = ProjectedDensityOfStates(fermi.*EV2HARTREE, pdos_raw[1,:].*EV2HARTREE, pdos_raw[2:end,:])
        end
    else
        # Don't return a vector filled with a bunch of undefined objects...
        pdos = ProjectedDensityOfStates[]
    end
    return (tdos, pdos)
end

readDOSCAR(file) =  open(readDOSCAR, convert_vasp_path(file, "DOSCAR"))
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
        energies .* EV2HARTREE,
        projband,
        phase
    )
end

readPROCAR(file) = open(readPROCAR, convert_vasp_path(file, "PROCAR"))
readPROCAR() = open(readPROCAR, "PROCAR")

"""
    get_fermi(file) -> NamedTuple{(:fermi, :alphabeta), NTuple{2,Float64}}

Reads a VASP OUTCAR file and returns the Fermi energy and alpha+beta value.
"""
function get_fermi(io::IO)
    readuntil(io, "E-fermi :")
    ln = split(readline(io))
    fermi = parse.(Float64, ln[1]) .* EV2HARTREE
    alphabeta = parse.(Float64, strip(ln[end],':')) .* EV2HARTREE
    return (fermi = fermi, alphabeta = alphabeta)
end

get_fermi(file) = open(get_fermi, convert_vasp_path(file, "OUTCAR"))

get_fermi() = open(get_fermi, "OUTCAR")
