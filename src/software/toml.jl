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

TOML.print(io::IO, x::AbstractAtomList; kwargs...) = TOML.print(io, to_toml_data(x); kwargs...)
TOML.print(io::IO, x::Crystal; kwargs...) = TOML.print(io, to_toml_data(x); kwargs...)
