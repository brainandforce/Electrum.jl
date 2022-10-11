#! /usr/bin/env julia

using Documenter
using Xtal

makedocs(
    sitename = "Xtal",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", nothing) == true)),
    modules = [Xtal],
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
    repo = "github.com/brainandforce/Xtal.jl"
)
