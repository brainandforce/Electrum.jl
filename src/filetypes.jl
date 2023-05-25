"""
    readXYZ(io::IO) -> AtomList{3}

Reads an XYZ file into an `AtomList{3}`.
"""
function readXYZ(io::IO)
    # Loop through each line in the file
    itr = eachline(io)
    # Get the number of atoms from the first line
    n_at = parse(Int, iterate(itr)[1])
    v = Vector{CartesianAtomPosition{3}}(undef, n_at)
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
        v[n] = CartesianAtomPosition{3}(name, pos)
    end
    return AtomList(v)
end

readXYZ(filename) = open(readXYZ, filename)

"""
    writeXYZ(io::IO, data::AbstractVector{<:AbstractAtomPosition})
    writeXYZ(io::IO, data::AbstractAtomList)
    writeXYZ(io::IO, data::AbstractCrystal)

Write an XYZ file based on a set of atomic coordinates.
"""
function writeXYZ(io::IO, data::AbstractVector{<:AbstractAtomPosition})
    # Write the number of atoms
    println(io, length(data))
    # Include the comment line
    println(io, "File written by Electrum.jl")
    # Write lines for all atoms
    for atom in data
        println(io, atomname(atom), join([lpad(@sprintf("%f", n), 11) for n in coord(atom)]))
    end
    return nothing
end

writeXYZ(io::IO, data::AbstractAtomList) = writeXYZ(io, coord(cartesian(data)))
writeXYZ(io::IO, data::AbstractCrystal) = writeXYZ(io, atoms(data))

function writeXYZ(filename, data) 
    open(filename; write=true) do io
        writeXYZ(io, data)
    end
    return nothing
end
