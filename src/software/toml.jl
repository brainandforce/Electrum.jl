function TOML.print(io::IO, b::AbstractBasis)
    println(io, "space = \"re", "cipro"^(DataSpace(b) isa ByReciprocalSpace), "al\"")
    println(io, "vectors = [", )
    for v in b
        print(io, " "^4, "[")
        join(io, v, ", ")
        println(io, "],")
    end
    println(io, "]")
end

function TOML.print(io::IO, p::PeriodicAtomList)
    println(io, "[basis]")
    TOML.print(io, basis(p))
    for atom in p
        println(io, "\n[[atoms]]")
        println(io, "name = ", name(atom))
        println(io, "number = ", atomic_number(atom))
        println(io, "coordinate = ", displacement(atom))
    end
end

function TOML.print(io::IO, xtal::Crystal)
    # Transform goes in the top-level table
    println(io, "transform = [")
    for v in eachcol(xtal.transform)
        print(io, " "^4, "[")
        join(io, v, ", ")
        println(io, "],")
    end
    println(io, "]\n")
    TOML.print(io, xtal.atoms)
    println(io, "\n[space_group]")
    println(io, "number = ", xtal.sgno)
    println(io, "origin = ", xtal.sgorig)
end
