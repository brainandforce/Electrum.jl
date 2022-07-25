using Documenter
using Xtal

makedocs(
    sitename = "Xtal",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", nothing) == true)),
    modules = [Xtal],
    pages = [
        "Home" => "index.md",
        "Types" => "types.md",
        "API" => Any[
            "Lattices" => "api/lattices.md"
            "Atoms" => "api/atoms.md"
            "Crystals" => "api/crystals.md"
            "Crystal data" => "api/data.md"
            "File types" => "api/filetypes.md"
        ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
