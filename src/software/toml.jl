"""
    Electrum.to_toml_data(x) -> Dict{String}

Converts an Electrum data type to a dictionary, which `TOML.print()` supports for direct output.
"""
function to_toml_data end

function to_toml_data(b::AbstractBasis)
    return Dict{String,Any}(
        "realspace" => !(DataSpace(b) isa ByReciprocalSpace),
        "dimension" => size(b,1),
        "vectors" => b.vectors
    )
end

to_toml_data(a::NamedAtom) = Dict{String,Any}("name" => name(a), "number" => atomic_number(a))

function to_toml_data(p::CartesianAtomPosition)
    data = to_toml_data(NamedAtom(p))
    data["cartesian"] = displacement(p)
    return data
end

function to_toml_data(p::FractionalAtomPosition)
    data = to_toml_data(NamedAtom(p))
    data["fractional"] = displacement(p)
    return data
end

function to_toml_data(p::FractionalAtomPosition, b::AbstractBasis)
    data = to_toml_data(p)
    data["cartesian"] = RealBasis(b) * displacement(p)
    return data
end

function to_toml_data(l::AtomList)
    return Dict{String,Any}(to_toml_data.(l))
end

function to_toml_data(l::PeriodicAtomList)
    return Dict{String,Any}("basis" => to_toml_data(basis(l)), "atoms" => to_toml_data.(l))
end

function to_toml_data(x::Crystal)
    data = to_toml_data(x.atoms)
    data["space_group"] = Dict{String,Any}("number" => x.sgno, "origin" => x.sgorig)
    data["transform"] = eachcol(x.transform)
    return data
end

writeTOML(io::IO, x::AbstractAtomList; kwargs...) = TOML.print(io, to_toml_data(x); kwargs...)
writeTOML(io::IO, x::Crystal; kwargs...) = TOML.print(io, to_toml_data(x); kwargs...)

"""
    writeTOML(file, x; sorted=false, by=identity)

Writes `x`, an Electrum data type, to a file. Currently, `x` may be either an `AbstractAtomList` or
a `Crystal`.
"""
function writeTOML(file, x; kwargs...) = open(io -> writeTOML(io, x; kwargs...), file)

export writeTOML
