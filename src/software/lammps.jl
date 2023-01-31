"""
    write_lammps_data(io::IO, list::PeriodicAtomList, [transform]; dummy::Bool = false)

Writes crystal information to a LAMMPS data format that can be used to define a simulation box
for running a molecular dynamics simulation. 
    
If `transform` is supplied, the list of atoms will be converted to a supercell with the associated
transformation (either a matrix, vector, or scalar). The `dummy` keyword determines whether dummy
atoms are included in the output (`false` by default).

This function currently only works for 3D systems.
"""
function write_lammps_data(io::IO, list::PeriodicAtomList; dummy::Bool=false) 
    println(
        io, 
        "# Written by Electrum.jl (https://github.com/brainandforce/Electrum.jl)"
    )
    # Get the number of atoms; write the corresponding line
    natoms = length(list) - length(filter(isdummy, list)) * !dummy
    println(io, natoms, " atoms")
    # Get the number of atom types
    println(io, natomtypes(list; dummy), " atom types")
    # Now print the new basis vectors
    println(io, "# Basis vector lengths:")
    println(io, @sprintf("0.000000    %f", basis(list)[1,1]), "     xlo xhi")
    println(io, @sprintf("0.000000    %f", basis(list)[2,2]), "     ylo yhi")
    println(io, @sprintf("0.000000    %f", basis(list)[3,3]), "     zlo zhi")
    # Add in tilt factors if needed
    if !isdiag(basis(list))
        println(io, "# Tilt factors:")
        @printf(io, "%f %f %f xy xz yz\n", basis(list)[1,2], basis(list)[1,3], basis(list)[2,3])
    end
    # Get the different atom types; strip dummy atoms if needed
    atl = (dummy ? list : filter(!isdummy, list))
    # Enumerate the atoms - by name, not atomic number
    names = name.(atomtypes(atl; dummy))
    # Print the atoms section
    println(io, "\nAtoms  # atomic\n")
    # Loop through the atoms
    for (n,atom) in enumerate(atl)
        @printf(
            io, "%i  %i  %f  %f  %f\n",
            # Atom number, enumerated atom type, coordinates (currently only 3D)...
            n, findfirst(isequal(name(atom)), names), (basis(atl) * atom.pos)...
        )
    end
end

function write_lammps_data(io::IO, list::PeriodicAtomList, transform; kwargs...)
    write_lammps_data(io, supercell(list, transform); kwargs...)
end

"""
    write_lammps_data(io::IO, xtal::AbstractCrystal, [transform]; dummy::Bool = false)

Writes crystal information to a LAMMPS data format that can be used to define a simulation box
for running a molecular dynamics simulation. 
    
The list of atoms that is written is given by converting `xtal` to a `PeriodicAtomList`, which uses
the supplied transformation matrix to generate all atomic positions. If `transform` is supplied,
the transformation will be applied to the `PeriodicAtomList` - it does not replace the transform
provided with `xtal`.

This function currently only works for 3D systems.
"""
function write_lammps_data(io::IO, xtal::AbstractCrystal;  kwargs...)
    write_lammps_data(io, PeriodicAtomList(xtal); kwargs...)
end

function write_lammps_data(io::IO, xtal::AbstractCrystal, transform; kwargs...)
    write_lammps_data(io, PeriodicAtomList(xtal), transform; kwargs...)
end

"""
    write_lammps_data(filename::AbstractString, data, [transform]; dummy::Bool = false)

Writes a LAMMPS data file to the path given by `filename`.
"""
function write_lammps_data(filename::AbstractString, args...; kwargs...)
    open(filename, write=true) do io
        write_lammps_data(io, args...; kwargs...)
    end
end
