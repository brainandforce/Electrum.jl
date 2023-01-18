#! /usr/bin/env julia

using Documenter
using Electrum

makedocs(
    sitename = "Electrum.jl",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", nothing) == true)),
    modules = [Electrum],
    pages = [
        "Home" => "index.md",
        "Types" => "types.md",
        "File formats" => "filetypes.md",
        "API" => Any[
            "Lattices" => "api/lattices.md"
            "Atoms" => "api/atoms.md"
            "Crystals" => "api/crystals.md"
            "Crystal data" => "api/data.md"
            "File formats" => "api/filetypes.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/brainandforce/Electrum.jl"
)
