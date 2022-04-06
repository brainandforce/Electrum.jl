"""
    readXYZ(io::IO) -> Vector{AtomPosition{3}}

Reads an XYZ file into a `Vector{AtomPosition{3}}`.
"""
function readXYZ(io::IO)
    # Loop through each line in the file
    itr = eachline(io)
    # Get the number of atoms from the first line
    n_at = parse(Int, iterate(itr)[1])
    v = Vector{AtomPosition{3}}(undef, n_at)
    # Skip the second line
    iterate(itr)
    # Loop through the rest
    for (n, ln) in enumerate(itr)
        # If the line is empty, skip it
        isempty(ln) && continue
        # Split up the current line
        entries = split(ln)
        # Get the name (first entry)
        name = entries[1]
        # Get the position (next three entries)
        pos = parse.(Float64, entries[2:4])
        # Add to vector
        v[n] = AtomPosition{3}(name, pos)
    end
    return v
end

readXYZ(filename::AbstractString) = open(readXYZ, filename)

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
        return BasisVectors{3}(hcat(vecs...))
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
            apos[n] = AtomPosition{3}(num, basis\pos)
        end
        return AtomList{3}(basis, apos)
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
    prim = zeros(BasisVectors{3})
    conv = zeros(BasisVectors{3})
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
            @debug "Line contents: \"" * ln * "\""
            # Loop until the end of the block
            while !contains(ln, "END_BLOCK")
                # Check for a datagrid entry
                if contains(ln, "DATAGRID_3D")
                    # Get the name of the datagrid
                    @debug "Line contents: \"" * ln * "\""
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
                    data[name] = RealSpaceDataGrid(latt, orig, grid)
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
    if iszero(conv)
        conv = BasisVectors{3}(REDUCTION_MATRIX_3D[ctr] \ matrix(prim))
    end
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

"""
    writeXSF(io, xtal::Crystal{D})

Writes the crystal component of an XCrysDen XSF file.
"""
function writeXSF(io::IO, xtal::Crystal{D}) where D
    println(io, "# Written by Xtal.jl")
    # Dimension information
    println(io, "DIM_GROUP")
    println(io, D, "  1")
    # Primitive cell vectors
    println(io, "PRIMVEC")
    for (n, x) in enumerate(prim(xtal))
        @printf(io, "%20.14f  ", x)
        n % D == 0 && println(io)
    end
    # Coordinates of the atoms in the primitive cell
    println(io, "PRIMCOORD")
    # Throw warning if the space group number is greater than 1
    # This probably means the generating set of atoms is incomplete
    xtal.sgno > 1 && @warn "Symmetrization of the unit cell has not been implemented yet."
    # Print the number of atoms in the structure
    println(io, lpad(string(natom(xtal.gen)), 6), "    1")
    # Print all of the atoms
    for a in cartesian(xtal.gen)
        @printf(io, "%6i", atomicno(a))
        for x in coord(a)
            @printf(io, "%20.14f  ", x)
        end
        println(io)
    end
end

"""
    writeXSF(io::IO, key, data::RealSpaceDataGrid{D,T}; periodic=true)

Writes the crystal component of an XCrysDen XSF file. By default, automatic wrapping of the 
datagrid occurs (values are repeated at the end of each dimension).
"""
function writeXSF(io::IO, key, data::RealSpaceDataGrid{D,T}; periodic=true) where {D,T}
    println(io, "    DATAGRID_", D, "D_", key)
    # Print the grid size (+1)
    println(io, " "^8, join([rpad(string(n + 1), 8) for n in gridsize(data)]))
    # Print the grid offset
    print(io, " "^4)
    for x in data.orig
        @printf(io, "%20.14f", x)
    end
    print(io, "\n" * " "^4)
    # Print the basis vectors for the grid
    for (n, x) in enumerate(basis(data))
        @printf(io, "%20.14f", x)
        n % D == 0 && print(io, "\n" * " "^4)
    end
    # This should generate a view that wraps around for a "general grid"
    if periodic
        g = view(grid(data), (1 .+ (0:n) .% n for n in gridsize(data))...)
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

function writeXSF(
    io::IO,
    data::Dict{K, RealSpaceDataGrid{D,T}};
    periodic=true
) where {K,D,T}
    println(io, "BEGIN_BLOCK_DATAGRID_", D, "D")
    println(io, "Written by Xtal.jl")
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

function writeXSF(filename::AbstractString, xtal::Crystal)
    open(filename, write=true) do io
        writeXSF(io, xtal)
    end
end

function writeXSF(filename::AbstractString, xtaldata::CrystalWithDatasets; periodic=true)
    open(filename, write=true) do io
        writeXSF(io, xtaldata, periodic=periodic)
    end
end

"""
    readCPcoeff(io::IO, Lmax::Val{L}) -> SphericalComponents{L}

Reads in the spherical harmonic projection coefficients from a CPpackage2 calculation.

By default, CPpackage2 gives the coefficients for spherical harmonics up to a maximum l value of 6.
"""
function readCPcoeff(io::IO, Lmax::Val{L}=Val(6)) where L
    # All the data should be in a Vector{Float64} with this one line
    data = parse.(Float64, [v[2] for v in split.(readlines(io))])
    @debug "$(length(data)) lines in file"
    # Number of spherical harmonic coefficients
    ncoeff = (L + 1)^2
    natom = div(length(data), ncoeff)
    @debug "ncoeff = $ncoeff, natom = $natom"
    return [SphericalComponents{L}(data[(n - 1)*ncoeff .+ (1:ncoeff)]) for n in 1:natom]
end

readCPcoeff(filename::AbstractString) = open(readCPcoeff, filename)

"""
    readCPgeo(io::IO)

Reads the basis vectors used for a CPpackage2 calculation.
"""
function readCPgeo(io::IO)

end

# Include code for processing specific types of files
include("abinit.jl")
include("vasp.jl")
