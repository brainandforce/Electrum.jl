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
    # Function for getting lattice basis vectors
    function getlattice!(itr)
        vecs = [parse.(Float64, split(iterate(itr)[1])) for n in 1:3]
        @debug string("Found vectors:\n", vecs)
        return hcat(vecs...)
    end
    # Function for getting 3D lists of atoms
    function getatoms!(itr, basis, natom)
        # Vector of AtomPositions
        apos = Vector{AtomPosition{3}}(undef, natom)
        for n in 1:natom
            # Get next line
            entries = split(iterate(itr)[1])
            # Get atomic number (first entry)
            num = parse(Int, entries[1])
            # Get atomic position (next 3)
            pos = SVector{3,Float64}(parse.(Float64, entries[s]) for s in 2:4)
            # assign atom position struct
            apos[n] = AtomPosition{3}(num, pos)
        end
        return AtomList{3}(basis, apos)
    end
    # Function to get grid data
    # By default, 
    function getgrid!(itr, sz; trim=true)
        # Store everything in this vector (to be reshaped)
        v = Vector{Float64}(undef, prod(sz))
        # Entry counter
        entry = 0
        # Line to be processed
        ln = ""
        while entry < prod(sz) || !contains(ln, "END_DATAGRID")
            # Get the next line
            ln = iterate(itr)[1]
            # Get the data entries
            dt = filter(!isnothing, [tryparse(Float64, s) for s in split(ln)])
            # Append the entries to the vector
            v[entry .+ (1:length(dt))] = dt
            # Move the entry counter over
            entry = entry + length(dt)
        end
        # Reshape and trim the edges of the grid, if needed
        return reshape(v, sz[1], sz[2], sz[3])[1:end - trim, 1:end - trim, 1:end - trim]
    end
    # Manually iterate through the file
    # This allows lines to be skipped
    iter = eachline(io)
    # Counter for debugging
    count = 0
    # Preallocated variables
    # Basis vectors
    prim = zeros(MMatrix{3,3,Float64})
    conv = zeros(MMatrix{3,3,Float64})
    # Atomic positions
    local atom_list::AtomList{3}
    data = Dict{String,RealSpaceDataGrid{3,Float64}}()
    # There is a break statement to get of out of this loop
    while true
        # Get the line contents
        ln = let i = iterate(iter)
            # Iteration returns nothing once it's out of lines, so break from loop
            isnothing(i) ? break : i[1]
        end
        count += 1
        # Check for empty lines or comment lines; advance if present
        # Comments are only allowed between blocks in a valid XSF file
        # No need to check within the loops below
        (isempty(ln) || startswith(ln, "#") || all(isspace, ln)) && continue
        # The block is supposed to be given as "BEGIN_BLOCK_DATAGRID_3D"
        # But bin2xsf seems to write "BEGIN_BLOCK_DATAGRID3D" (missing underscore)
        if contains(ln, "BEGIN_BLOCK_DATAGRID") && contains(ln, "3D")
            # Skip this line
            ln = iterate(iter)[1]
            # Loop until the end of the block
            while !contains(ln, "END_BLOCK")
                ln = iterate(iter)[1]
                # Check for a datagrid entry
                if contains(ln, "DATAGRID_3D")
                    # Get the name of the datagrid
                    name = split(ln, "DATAGRID_3D_")[2]
                    # Get the dimensions of the datagrid
                    ln = iterate(iter)[1]
                    dims = parse.(Int, split(ln))
                    grid = Array{Float64,3}(undef, dims...)
                    # Get the shift off the origin
                    ln = iterate(iter)[1]
                    orig = parse.(Float64, split(ln))
                    # Get the lattice
                    latt = getlattice!(iter)
                    # Get the datagrid
                    # TODO: Determine if trimming is needed based on periodicity
                    # Periodic structures should get trimmed
                    grid = getgrid!(iter, dims, trim=true)
                    # Export the data in a RealSpaceDataGrid
                    data[name] = RealSpaceDataGrid{3,Float64}(latt, orig, grid)
                end
            end
        # Get the primitive cell vectors
        elseif contains(ln, "PRIMVEC")
            prim = getlattice!(iter)
        # Get the conventional cell vectors (if present)
        elseif contains(ln, "CONVVEC")
            conv = getlattice!(iter)
        # Get primitive cell atomic coordinates
        elseif contains(ln, "PRIMCOORD")
            # Get the number of atoms
            ln = iterate(iter)[1]
            natom = parse(Int, split(ln)[1])
            # The next lines have all of the atoms
            atom_list = getatoms!(iter, prim, natom)
        # Get conventional cell atomic coordinates
        # But do we need them?
        elseif contains(ln, "CONVCOORD") && isempty(atom_list)
            # TODO: complete this, but it's low priority
            # The PRIMCOORD block gets everything we need right now
        elseif contains(ln, "ATOMS") && isempty(atom_list)
            # TODO: complete this, but it's low priority
            # The PRIMCOORD block gets everything we need right now
        end
    end
    # If the conventional cell hasn't been defined, generate it
    iszero(conv) && (conv = prim * inv(REDUCTION_MATRIX_3D[ctr]))
    # Generate the real space lattice
    latt = RealLattice{3}(prim, conv)
    xtal = Crystal{3}(latt, spgrp, origin, atom_list, atom_list)
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

# Include code for processing specific types of files
include("abinit.jl")
include("vasp.jl")