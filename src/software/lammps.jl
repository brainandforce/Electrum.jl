"""
    write_lammps_data(
        io::IO,
        xtal::AbstractCrystal{D};
        supercell = diagm(ones(Int, D))::AbstractVector{<:Integer};
        dummy::Bool = false 
    )

Writes crystal information to a LAMMPS data format that can be used to define a simulation box
for running a molecular dynamics simulation.

This function currently only works for 3D systems.
"""
function write_lammps_data(
    io::IO,
    list::AtomList{D},
    supercell::AbstractMatrix{<:Integer};
    dummy::Bool=false
) where D
    println(io, "# LAMMPS data file written by Xtal.jl (https://github.com/brainandforce/Xtal.jl)")
    # Get the number of atoms; write the corresponding line
    natoms = length(list) * convert(Int, abs(det(supercell)))
    println(io, natoms, " atoms")
    # Get the number of atom types
    println(io, natomtypes(list), " atom types")
    # Convert to upper triangular form
    bnew = triangularize(basis(list), supercell)
    # Now print the new basis vectors
    println(io, @sprintf("0.000000    %f", bnew[1,1]), "     xlo xhi")
    println(io, @sprintf("0.000000    %f", bnew[2,2]), "     ylo yhi")
    println(io, @sprintf("0.000000    %f", bnew[3,3]), "     zlo zhi")
    # Add in tilt factors if needed
    if !isdiag(matrix(bnew))
        println(io, @sprintf("%f %f %f", bnew[1,2], bnew[1,3], bnew[2,3]), " xy xz yz")
    end
    # Get the different atom types
    # Strip dummy atoms if needed
    atl = (dummy ? list : remove_dummies(list))
    # Enumerate the atoms - by name, not atomic number
    names = atomnames(atl)
    # Loop through the atoms
    for (n,atom) in enumerate(atl)
        @printf(
            io, "%i  %i  %f  %f  %f",
            # Atom number, enumerated atom type, coordinates (currently only 3D)...
            n, findfirst(isequal(atomname(atom)), names), (basis(atl) * atom.pos)...
        )
    end
end

# Same thing, but with an optional supercell vector
# (This is automatically converted to a diagonal matrix)
write_lammps_data(
    io::IO, 
    list::AtomList{D},
    supercell::AbstractVector{<:Integer} = ones(Int, D);
    kwargs...
) where D = write_lammps_data(io, list, diagm(supercell); kwargs...)

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
    xtal::Union{AbstractCrystal{D},AtomList{D}}, # this is weird...
    supercell = ones(Int, D);
    kwargs...
) where D
    open(filename, write=true) do io
        write_lammps_data(io, xtal, supercell; kwargs...)
    end
end
