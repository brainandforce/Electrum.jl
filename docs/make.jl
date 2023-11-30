#! /usr/bin/env julia

using Documenter
using LinearAlgebra, FFTW
using Electrum

is_ci_env = (get(ENV, "CI", nothing) == true)
@info "is_ci_env == $is_ci_env"

makedocs(
    sitename = "Electrum.jl",
    format = Documenter.HTML(prettyurls = is_ci_env),
    modules = [Electrum],
    checkdocs = :exports,
    warnonly = is_ci_env,
    pages = [
        "Home" => "index.md",
        "Lattices" => "lattices.md",
        "Atoms and crystals" => "atoms.md",
        "Data grids" => "grids.md",
        "File formats" => "filetypes.md",
        "API" => Any[
            "Geometry" => "api/geometry.md"
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
