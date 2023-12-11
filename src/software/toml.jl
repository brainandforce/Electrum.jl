"""
    Electrum.toml_convert(x) -> Dict{String}

Converts an Electrum data type to a dictionary, which `TOML.print()` supports for direct output.
"""
function toml_convert end

function toml_convert(b::LatticeBasis)
    return Dict{String,Any}(
        "realspace" => !(BySpace(b) === ByReciprocalSpace()),
        "dimension" => size(b, 1),
        "vectors" => b.vectors
    )
end

toml_convert(a::NamedAtom) = Dict{String,Any}("name" => name(a), "number" => atomic_number(a))

function toml_convert(p::CartesianAtomPosition)
    data = toml_convert(NamedAtom(p))
    data["cartesian"] = displacement(p)
    return data
end

function toml_convert(p::FractionalAtomPosition)
    data = toml_convert(NamedAtom(p))
    data["fractional"] = displacement(p)
    return data
end

function toml_convert(p::FractionalAtomPosition, b::LatticeBasis)
    data = toml_convert(p)
    data["cartesian"] = RealBasis(b) * displacement(p)
    return data
end

toml_convert(l::AtomList) = toml_convert.(l)

function toml_convert(l::PeriodicAtomList)
    return Dict{String,Any}("basis" => toml_convert(basis(l)), "atoms" => toml_convert.(l))
end

function toml_convert(x::Crystal{D}) where D
    data = toml_convert(x.atoms)
    data["space_group"] = Dict{String,Any}("number" => x.sgno, "origin" => x.sgorig)
    data["transform"] = SVector{D}(eachcol(x.transform))
    return data
end

writeTOML(io::IO, x::AbstractAtomList; kwargs...) = TOML.print(io, toml_convert(x); kwargs...)
writeTOML(io::IO, x::Crystal; kwargs...) = TOML.print(io, toml_convert(x); kwargs...)

"""
    writeTOML(file, x; sorted=false, by=identity)

Writes `x`, an Electrum data type, to a file. Currently, `x` may be either an `AbstractAtomList` or
a `Crystal`.

The `sorted` and `by` keywords are passed through from `TOML.print` and allow for sorting of the
keys in the file output.
"""
writeTOML(file, x; kwargs...) = open(io -> writeTOML(io, x; kwargs...), file; write=true)

export writeTOML
