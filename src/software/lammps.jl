"""
    read_lammps_data(io::IO; atoms::AbstractVector = NamedAtom[]) -> PeriodicAtomList{3}

Reads a LAMMPS data file containing atomic coordinates. Currently, this function only supports 3D
crystal data.
"""
function read_lammps_data(io::IO; atom_types::AbstractVector = NamedAtom[]) 
    natom = 0
    ntypat = length(atom_types)
    basis_M = zeros(MMatrix{3,3,Float64})
    while !eof(io)
        # Get rid of any comment characters
        ln = strip(first(split(readline(io), '#')))
        # Get the number of atoms
        if isempty(ln)
            continue
        elseif ln == "Atoms"
            break
        elseif occursin(" atoms", ln)
            natom = parse(Int, first(split(ln)))
        elseif occursin(" atom types", ln)
            ntypat = parse(Int, first(split(ln)))
        elseif occursin(r"xlo.+xhi", ln)
            basis_M[1,1] = parse(Float64, split(ln)[2]) - parse(Float64, split(ln)[1])
        elseif occursin(r"ylo.+yhi", ln)
            basis_M[2,2] = parse(Float64, split(ln)[2]) - parse(Float64, split(ln)[1])
        elseif occursin(r"zlo.+zhi", ln)
            basis_M[3,3] = parse(Float64, split(ln)[2]) - parse(Float64, split(ln)[1])
        # Match these in any order
        elseif all(occursin(t, ln) for t in ["xy", "xz", "yz"])
            tilt_tokens = filter(s -> (s in ["xy", "xz", "yz"]), split(ln))
            indices = [
                Tuple(parse(Int, c) for c in s)
                for s in replace.(tilt_tokens, "x" => "1", "y" => "2", "z" => "3")
            ]
            for (i,s) in zip(indices, split(ln)[1:3])
                basis_M[i...] = parse.(Float64, s)
            end
        end
    end
    atom_vec = Vector{FractionalAtomPosition{3}}(undef, natom)
    ind = 1
    while ind <= natom && !eof(io)
        # Get rid of any comment characters
        ln = strip(first(split(readline(io), '#')))
        isempty(ln) && continue
        # Convert the Cartesian coordinate to fractional
        coord_cart = SVector{3,Float64}(parse(Float64, s) for s in split(ln)[3:end])
        coord_frac = basis_M \ coord_cart
        atom_vec[ind] = if !isempty(atom_types)
            FractionalAtomPosition(atom_types[parse(Int, split(ln)[2])], coord_frac)
        else
            FractionalAtomPosition(0, coord_frac)
        end
        ind += 1
    end
    return PeriodicAtomList(RealBasis(basis_M), atom_vec)
end

read_lammps_data(filename; kwargs...) = open(io -> read_lammps_data(io; kwargs...), filename)

"""
    write_lammps_data(io::IO, list::PeriodicAtomList, [transform]; dummy::Bool = false)

Writes crystal information to a LAMMPS data format that can be used to define a simulation box for
running a molecular dynamics simulation. 
    
If `transform` is supplied, the list of atoms will be converted to a supercell with the associated
transformation (either a matrix, vector, or scalar). The `dummy` keyword determines whether dummy
atoms are included in the output (`false` by default).

This function currently only works for 3D systems.
"""
function write_lammps_data(io::IO, list::PeriodicAtomList{D}; dummy::Bool=false) where D
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
    isdiag(basis(list)) || @printf(
        io, "# Tilt factors:\n%f %f %f xy xz yz\n",
        basis(list)[1,2], basis(list)[1,3], basis(list)[2,3]
    )
    # Get the different atom types; strip dummy atoms if needed
    atl = (dummy ? list : filter(!isdummy, list))
    # Enumerate the atoms - by name, not atomic number
    names = name.(atomtypes(atl; dummy))
    # Print the atoms section
    println(io, "\nAtoms  # atomic\n")
    # Determine the correct padding amounts for the leftmost digits
    atomno_pad = length(digits(length(list))) + 2
    atomtype_pad = length(digits(natomtypes(list)))
    # Loop through the atoms
    for (n,atom) in enumerate(atl)
        print(io, rpad(n, atomno_pad))
        print(io, rpad(findfirst(isequal(name(atom)), names), atomtype_pad))
        for m in 1:D
            @printf(io, "% 16.10f", (basis(atl) * displacement(atom))[m])
        end
        println(io)
    end
end

function write_lammps_data(io::IO, list::PeriodicAtomList, transform; kwargs...)
    write_lammps_data(io, supercell(list, transform); kwargs...)
end

"""
    write_lammps_data(io::IO, xtal::AbstractCrystal, [transform]; dummy::Bool = false)

Writes crystal information to a LAMMPS data format that can be used to define a simulation box for
for running a molecular dynamics simulation. 
    
The list of atoms that is written is given by converting `xtal` to a `PeriodicAtomList`, which uses
the supplied transformation matrix to generate all atomic positions. If `transform` is supplied, the
transformation will be applied to the `PeriodicAtomList` - it does not replace the transform
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
    write_lammps_data(file, data, [transform]; dummy::Bool = false)

Writes a LAMMPS data file to the path given by `filename`.
"""
function write_lammps_data(filename, args...; kwargs...)
    open(filename, write=true) do io
        write_lammps_data(io, args...; kwargs...)
    end
end
