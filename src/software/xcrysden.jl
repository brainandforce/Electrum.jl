"""
    readXSF3D(
        io::IO;
        spgrp::Integer = 0,
        origin::AbstractVector{<:Real} = [0, 0, 0]
        ctr::Symbol = :P
    ) -> CrystalWithDatasets{3}

Reads in an XCrysDen XSF file from an input stream and returns a `CrystalWithDatasets{3}` with all
datasets that have been included within the file.

Space group and origin information are not supplied in XSF files, but they can be supplied using the
`spgrp` and `origin` keyword arguments. Centering information can be provided using the `ctr` 
argument, but is overridden by a space group assignment.
"""
function readXSF3D(
    io::IO;
    spgrp::Integer = 0,
    origin::AbstractVector{<:Real} = [0, 0, 0],
)
    # Function for getting lattice basis vectors
    function getlattice!(itr)
        latt = SVector{3}((SVector{3}(parse.(Float64, split(iterate(itr)[1]))) for _ in 1:3))
        @debug string("Found vectors:\n", (string(v) * "\n" for v in latt)...)
        return RealBasis(hcat(latt...) * ANG2BOHR)
    end
    # Function for getting 3D lists of atoms
    function getatoms!(itr, basis::LatticeBasis, natom::Integer)
        # Vector of AtomPositions
        apos = Vector{CartesianAtomPosition{3}}(undef, natom)
        for n in 1:natom
            # Get next line
            entries = split(iterate(itr)[1])
            # Get atomic number (first entry)
            num = parse(Int, entries[1])
            # Get atomic position (next 3)
            pos = SVector{3,Float64}(parse.(Float64, entries[s]) for s in 2:4)
            # Add this to the vector
            apos[n] = CartesianAtomPosition{3}(num, pos)
        end
        return PeriodicAtomList(basis, apos)
    end
    # Function to get grid data
    # By default, trim the edges of the grid (repeated points)
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
    prim = zeros(RealBasis{3})
    conv = zeros(RealBasis{3})
    # Atomic positions
    local atom_list::PeriodicAtomList{3}
    data = Dict{String,RealDataGrid{3,Float64}}()
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
            @debug "Line contents: \"" * ln * "\""
            # Loop until the end of the block
            while !contains(ln, "END_BLOCK")
                # Check for a datagrid entry
                if contains(ln, "DATAGRID_3D")
                    # Get the name of the datagrid
                    @debug "Line contents: \"" * ln * "\""
                    name = rstrip(split(ln, "DATAGRID_3D_")[2])
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
                    # Export the data in a RealDataGrid
                    data[name] = RealDataGrid(grid, latt, orig)
                end
                ln = iterate(iter)[1]
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
            @debug string("natom = ", natom)
            # The next lines have all of the atoms
            atom_list = getatoms!(iter, prim, natom)
        elseif contains(ln, "CONVCOORD") && isempty(atom_list)
            # TODO: complete this, but it's low priority
            # The PRIMCOORD block gets everything we need right now
        elseif contains(ln, "ATOMS") && isempty(atom_list)
            # TODO: complete this, but it's low priority
            # The PRIMCOORD block gets everything we need right now
        elseif contains(ln, "DIM-GROUP")
            # Check the line below and see if the structure is periodic
            ln = iterate(iter)[1]
            if first(strip(ln)) == '3' && spgrp == 0
                # If the structure is a crystal, set spgrp to 1
                # spgrp 0 implies a non-periodic structure
                spgrp = 1
            end
        elseif contains(ln, "CRYSTAL") && spgrp == 0
            # If the structure is a crystal, set spgrp to 1
            # spgrp 0 implies a non-periodic structure
            spgrp = 1
        end
    end
    # Generate the transform between primitive and conventional lattices, if needed
    transform = if iszero(conv) 
        SMatrix{3,3,Int}(LinearAlgebra.I)
    else
        round.(Int, conv.matrix / prim.matrix)
    end
    return CrystalWithDatasets{3,String,RealDataGrid{3,Float64}}(
        Crystal(atom_list, spgrp, origin, transform), data
    )
end

"""
    readXSF3D(
        file;
        spgrp::Integer = 0,
        origin::AbstractVector{<:Real} = [0, 0, 0]
        ctr::Symbol = :P
    ) -> CrystalWithDatasets{3}

Reads an XSF file at path `filename`.
"""
function readXSF3D(filename; kwargs...)
    open(filename) do io
        readXSF3D(io; kwargs...)
    end
end

const readXSF = readXSF3D

"""
    writeXSF(io, xtal::PeriodicAtomList{D})

Writes the crystal component of an XCrysDen XSF file.
"""
function writeXSF(io::IO, l::PeriodicAtomList{D}) where D
    println(io, "# Written by Electrum.jl")
    # Dimension information
    println(io, "DIM_GROUP")
    println(io, D, "  1")
    # Primitive cell vectors
    println(io, "PRIMVEC")
    for (n, x) in enumerate(basis(l))
        @printf(io, "%20.14f  ", x * BOHR2ANG)
        n % D == 0 && println(io)
    end
    # Coordinates of the atoms in the primitive cell
    println(io, "PRIMCOORD")
    # Print the number of atoms in the structure
    println(io, lpad(string(length(l)), 6), "    1")
    # Print all of the atoms
    for a in AtomList(l)
        @printf(io, "%6i", atomic_number(a))
        for x in displacement(a)
            @printf(io, "%20.14f  ", x)
        end
        println(io)
    end
end

writeXSF(io::IO, xtal::Crystal) = writeXSF(io, xtal.atoms)

"""
    writeXSF(io::IO, key, data::RealDataGrid{D}; periodic=true)

Writes the crystal component of an XCrysDen XSF file. By default, automatic wrapping of the datagrid
occurs (values are repeated at the end of each dimension).
"""
function writeXSF(io::IO, key, data::RealDataGrid{D}; periodic=true) where D
    println(io, "    DATAGRID_", D, "D_", key)
    # Print the grid size (+1)
    println(io, " "^8, join([rpad(string(n + 1), 8) for n in size(data)]))
    # Print the grid offset
    print(io, " "^4)
    for x in data.orig
        @printf(io, "%20.14f", x)
    end
    print(io, "\n" * " "^4)
    # Print the basis vectors for the grid
    for (n, x) in enumerate(matrix(basis(data)))
        @printf(io, "%20.14f", x)
        n % D == 0 && print(io, "\n" * " "^4)
    end
    # This should generate a view that wraps around for a "general grid"
    if periodic
        g = view(data.data, (1 .+ (0:n) .% n for n in size(data))...)
    else
        g = view(M, (1:n for n in size(M))...)
    end
    # Print the actual data in the grid
    for (n,x) in enumerate(g)
        @printf(io, "%20.14f", x)
        n % 4 == 0 && print(io, "\n" *  " "^4)
    end
    # End the datagrid
    length(g) % 4 == 0 || println(io)
    println(io, "    END_DATAGRID_", D, "D")
end

function writeXSF(io::IO, data::Dict{K,<:RealDataGrid{D}}; periodic=true) where {K,D}
    println(io, "BEGIN_BLOCK_DATAGRID_", D, "D")
    println(io, "Written by Electrum.jl")
    for (k,d) in data
        writeXSF(io::IO, k, d, periodic=periodic)
    end
    println(io, "END_BLOCK_DATAGRID_", D, "D")
end

"""
    writeXSF(io::IO, xtaldata::CrystalWithDatasets{D,K,V})

Writes a `CrystalWithDatasets` to an XCrysDen XSF file.
"""
function writeXSF(
    io::IO,
    xtaldata::CrystalWithDatasets;
    periodic=true
)
    writeXSF(io, Crystal(xtaldata))
    writeXSF(io, data(xtaldata), periodic=periodic)
end

function writeXSF(filename, data...; kwargs...)
    open(filename, write=true) do io
        writeXSF(io, data...; kwargs...)
    end
end
