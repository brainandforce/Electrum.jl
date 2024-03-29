# Electrum.jl

[![Documentation (stable)][docs-stable-img]][docs-stable-url]
[![Documentation (dev)][docs-dev-img]][docs-dev-url]
[![CI status][ci-status-img]][ci-status-url]
[![Codecov][codecov-img]][codecov-url]
[![Aqua.jl][aqua-img]][aqua-url]

A Julia package for working with crystal structures, associated data, and various file formats,
with the aim of making theory development for solid state chemistry and materials science easier
for everyone.

This package is written in pure Julia with minimal dependencies, so instead of pulling operations
from established packages in other languages (such as scipy) they are implemented in this package.
We hope that Electrum.jl can serve as a reference for both useful operations on crystal structures
and well-maintained scientific software.

## How to install

Electrum.jl is available in the Julia General package registry! You can add it to your environment
with the `add` command in package mode:
```
(@v1.6+) pkg> add https://github.com/brainandforce/Electrum.jl 
```
Or you can import the `Pkg` module:
```julia-repl
julia> Pkg.add("Electrum")
```
If you'd like to track a specific on this repository for development purposes, you can refer to the
repo URL and tag it with the desired branch name (here, I've used `release`):
```
(@v1.6+) pkg> add https://github.com/brainandforce/Electrum.jl#release
```
```julia-repl
julia> Pkg.add(url="https://github.com/brainandforce/Electrum.jl", rev="release")
```
The current development state is kept in the `main` branch, and the most recent stable version is
the head of the `release` branch. Specific releases for a minor version may be found by suffixing
the minor version with `/release`: for instance, `0.1/release` contains the latest release version
in the 0.1 series.

Electrum.jl is tested against the most recent LTS release (currently 1.6.7) and the current stable
release. Julia 1.5 and earlier are not supported.

## Current features

* Reading of common file formats:
     + abinit potential, density, and wavefunction outputs: tested on versions 7.10.5, 8.10.3,
    and 9.10.1 (current version as of last update)
     + LAMMPS atomic position data: tested on 3 November 2022 release
     + VASP POSCAR/CONTCAR, WAVECAR, DOSCAR, and PROCAR: tested on version 4.6
     + XCrysDen XSF: tested against [the official specification][xsf-spec-url]
     + XYZ files
     + CPpackage outputs: tested against the release version of [CPpackage3][cppackage-url]
* Writing of common file formats:
     + LAMMPS atomic position data
     + XCrysDen XSF
     + XYZ files
     + VASP POSCAR/CONTCAR
     + TOML files
* Operations on datagrids:
     + Standard arithmetic operations: addition, subtraction, multiplication, division, negation...
     + Broadcasting of any Julia operation with dot syntax
     + FFTs and inverse FFTs on real and reciprocal space data
     + Real space gradients

## Planned features

This project is just starting to get off the ground, but here's what we have planned:

  * Reading and writing of common file formats:
      + XTL files
      + CIF files
  * Manipulation of data grids associated with crystal structures
      + Real space grid reinterpolation
  * Other operations
      + Working with crystal symmetry
      + Converting k-point meshes to symmetry-reduced lists of weighted k-points
  * Integration with [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl)
     
...and more that we will decide in time! If you'd like to contibute, be sure to read the included
[contributing guidelines.](CONTRIBUTING.md)

## Packages which use Electrum.jl

If your package uses Electrum.jl, we'd like to know about it - open an issue or a pull request with
changes to the list of packages we know depend on Electrum.jl given below:
  - [DFT-raMO][dftramo-url] (not currently available in the Julia package registry)

[docs-stable-img]:  https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]:  https://brainandforce.github.io/Electrum.jl/stable
[docs-dev-img]:     https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:     https://brainandforce.github.io/Electrum.jl/dev
[ci-status-img]:    https://github.com/brainandforce/Electrum.jl/workflows/CI/badge.svg
[ci-status-url]:    https://github.com/brainandforce/Electrum.jl/actions
[aqua-img]:         https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]:         https://github.com/JuliaTesting/Aqua.jl
[codecov-img]:      https://codecov.io/gh/brainandforce/Electrum.jl/branch/main/graph/badge.svg
[codecov-url]:      https://codecov.io/gh/brainandforce/Electrum.jl/
[xsf-spec-url]:     http://www.xcrysden.org/doc/XSF.html
[cppackage-url]:    https://github.com/dcfredrickson/CPpackage3
[dftramo-url]:      https://github.com/xamberl/DFT-raMO
