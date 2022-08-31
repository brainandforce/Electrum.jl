"""
    write_lammps_data(
        io::IO,
        list::AtomList{D}
        sc::AbstractVecOrMat{<:Integer} = ones(Int,D);
        dummy::Bool = false
    )

Writes crystal information to a LAMMPS data format that can be used to define a simulation box
for running a molecular dynamics simulation. If provided, `sc` can be used to generate a supercell

This function currently only works correctly for 3D systems.
"""
function write_lammps_data(io::IO, list::AtomList{D}, dummy::Bool=false) where D
    println(io, "# LAMMPS data file written by Xtal.jl (https://github.com/brainandforce/Xtal.jl)")
    # Get the number of atoms; write the corresponding line
    natoms = length(list)
    println(io, natoms, " atoms")
    # Get the number of atom types
    println(io, natomtypes(list), " atom types")
    # Now print the new basis vectors
    println(io, @sprintf("0.000000    %f", bnew[1,1]), "     xlo xhi")
    println(io, @sprintf("0.000000    %f", bnew[2,2]), "     ylo yhi")
    println(io, @sprintf("0.000000    %f", bnew[3,3]), "     zlo zhi")
    # Add in tilt factors if needed
    if !isdiag(matrix(bnew))
        println(io, "# Tilt factors:")
        @printf(io, "%f %f %f xy xz yz", bnew[1,2], bnew[1,3], bnew[2,3])
    end
    # Get the different atom types; strip dummy atoms if needed
    atl = (dummy ? list : remove_dummies(list))
    # Enumerate the atoms - by name, not atomic number
    names = atomnames(atl)
    # Print the atoms section
    println(io, "\nAtoms  # atomic\n")
    # Loop through the atoms
    for (n,atom) in enumerate(atl)
        @printf(
            io, "%i  %i  %f  %f  %f\n",
            # Atom number, enumerated atom type, coordinates (currently only 3D)...
            n, findfirst(isequal(atomname(atom)), names), (basis(atl) * atom.pos)...
        )
    end
end

"""
    write_lammps_data(
        io::IO,
        xtal::AbstractCrystal{D};
        sc = diagm(ones(Int, D))::AbstractVector{<:Integer};
        dummy::Bool = false 
    )

Writes crystal information to a LAMMPS data format that can be used to define a simulation box
for running a molecular dynamics simulation.

This function currently only works for 3D systems.
"""
function write_lammps_data(
    io::IO,
    list::AtomList{D},
    sc::AbstractVecOrMat{<:Integer};
    kwargs...
) where D
    write_lammps_data(io, supercell(list, sc); kwargs...)
end

# Same thing, but with an AbstractCrystal
write_lammps_data(
    io::IO, 
    xtal::AbstractCrystal{D},
    supercell;
    kwargs...
) where D = write_lammps_data(io, xtal.gen, supercell; kwargs...)

# Write to a filename
function write_lammps_data(
    filename::AbstractString,
    xtal,
    supercell = ones(Int, D);
    kwargs...
) where D
    open(filename, write=true) do io
        write_lammps_data(io, xtal, supercell; kwargs...)
    end
end
