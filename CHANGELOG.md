# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

(Well, to be honest, we do the best we can - we're a bunch of relatively new developers and we're
trying our best to learn good practices as we build the tools we need for our research. Changes from
the 0.1 series aren't properly documented, but we're moving towards using strict changelogs and 
semantic versioning from 0.2 onwards.)

# Changes from the 0.1 series

## New features

* A documentation suite has been added using
[Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).
* A testing suite has been added using [Test.jl](https://docs.julialang.org/en/v1/stdlib/Test/).
* Citation is facilitated with [`CITATION.cff`](CITATION.cff).

## Breaking changes

### Basis vectors and lattices

* `AbstractLattice` and its subtypes (`RealLattice` and `ReciprocalLattice`) have been removed.
* `BasisVectors` has been removed in favor of `RealBasis` and `ReciprocalBasis`, which are subtypes
of `AbstractBasis`. These have the advantage of explicitly declaring whether the vectors describe
real or reciprocal space.

### Crystal data
* The `pos` field of `Crystal`, which was used to store explicit atomic positions, has been 
removed, as have all methods that depend on this field.

### Pseudopotentials
* All potential and pseudopotential functionality has been removed from this package and will be 
placed in a new package.

### Renaming of structs, functions, and constants
* `cell_angle_cos()`, `cell_angle_rad()`, and `cell_angle_deg()` have been simplified to 
`angles_cos()`, `angles_rad()`, and `angles_deg()`.
* 
