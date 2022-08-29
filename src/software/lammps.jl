"""
    write_lammps_data(
        io::IO,
        xtal::Crystal{D};
        supercell = diagm(ones(Int, D))::AbstractVector{<:Integer}
    )

Writes crystal information to a LAMMPS data format that can be used to define a simulation box
for running a molecular dynamics simulation.

This function currently only works for 3D systems.
"""
function write_lammps_data(
    io::IO,
    xtal::AbstractCrystal{D},
    supercell::AbstractMatrix{<:Integer};
    dummy::Bool=false
) where D
    println(io, "# LAMMPS data file written by Xtal.jl (https://github.com/brainandforce/Xtal.jl)")
    # Get the number of atoms; write the corresponding line
    natoms = natom_cell(xtal) * convert(Int, abs(det(supercell)))
    println(io, natoms, " atoms")
    # Get the number of atom types
    println(io, natomtypes(xtal), " atom types")
    # Convert to upper triangular form
    bnew = triangularize(basis(xtal, primitive=true), supercell)
    # Now print the new basis vectors
    println(io, @sprintf("0.000000    %f", bnew[1,1]), "     xlo xhi")
    println(io, @sprintf("0.000000    %f", bnew[2,2]), "     ylo yhi")
    println(io, @sprintf("0.000000    %f", bnew[3,3]), "     zlo zhi")
    # Add in tilt factors if needed
    if !isdiag(matrix(bnew))
        println(io, @sprintf("%f %f %f", bnew[1,2], bnew[1,3], bnew[2,3]), " xy xz yz")
    end
    # Get the different atom types
    # First enumerate the atoms - by name, not atomic number
    atl = (dummy ? xtal.gen : remove_dummies(xtal.gen))
    names = atomnames(atl)
    for (n,atom) in enumerate(atl)
        println(
            io, n, " ", findfirst(isequal(atomname(atom)), names),
            @sprintf("%f  %f  %f", basis(atl) * atom.pos)
        )
    end 
    # Print the atom list
    println(io, "\nAtoms  # atomic\n")
    for (n,atom) in enumerate(xtal.gen)
        println(io, n, " ", )
    end
end

write_lammps_data(
    io::IO, 
    xtal::AbstractCrystal{D},
    supercell::AbstractVector{<:Integer} = ones(Int, D)
) where D = write_lammps_data(io, xtal, diagm(supercell))

function write_lammps_data(filename::AbstractString, xtal::Crystal{D}; kwargs...) where D
    open(filename, write=true) do io
        write_lammps_data(io, xtal, kwargs...)
    end
end
