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

### Crystals
* The `Crystal` datatype has been completely refactored. Instead of consisting of a `RealLattice`,
the basis vectors of the generating atomic positions are used, and the `transform` field contains a
matrix describing how the basis vectors may be transformed to a conventional cell.
* The `pos` field of `Crystal`, containing a list of explicitly generated atomic positions, has
been removed. This was originally intended to be included to facilitate plotting, but this will
instead be added to a different package.

### Crystal data
* `RealSpaceDataGrid` now uses zero-based indexing so that data placed at `[0, 0, 0]` in fractional
coordinates can be accessed with indices `[0, 0, 0]`.
* `HKLData` now stores bounds information implicitly by using the natural FFT frequency binning as
indices.
* Both types now support zero-based linear indexing.
* `ReciprocalWavefunction` now stores its data in 3D matrices of size `spins * kpts * bands`, and
may now be indexed that way. Indexing a `ReciprocalWavefunction` returns a `NamedTuple` with keys
`:coeffs`, `:energies`, and `:occupancies`.

### Pseudopotentials
* All potential and pseudopotential functionality has been removed from this package and will be 
placed in a new package.

### Renaming of structs, functions, and constants
* `gridsize()` has been removed in favor of `Base.size()`.
* `lattice2d()` and `lattice3d()` have been removed
* `cell_angle_cos()`, `cell_angle_rad()`, and `cell_angle_deg()` have been simplified to 
`angles_cos()`, `angles_rad()`, and `angles_deg()`.
* `SphericalComponents` is now `SphericalHarmonic`.
