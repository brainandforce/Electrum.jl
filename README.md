# Xtal.jl

[![Documentation (stable)][docs-stable-img]][docs-stable-url]
[![Documentation (dev)][docs-dev-img]][docs-dev-url]
[![Aqua QA][aqua-img]][aqua-url]

A Julia package for working with crystal structures, associated data, and various file formats,
with the aim of making theory development for solid state chemistry and materials science easier
for everyone.

This package is written in pure Julia with minimal dependencies, so instead of pulling operations
from established packages in other languages (such as scipy) they are implemented in this package.
We hope that Xtal.jl can serve as a reference for both useful operations on crystal structures and
well-maintained scientific software.

## How to install

Xtal.jl is not in the Julia package registry yet, so you'll need to install it like this:
```
(@v1.8) pkg> add https://github.com/brainandforce/Xtal.jl 
```
You can access package mode by typing `]` at the REPL. Alternatively, you can work with the `Pkg`
module:
```julia-repl
julia> Pkg.add(url="https://github.com/brainandforce/Xtal.jl")
```
If you'd like to track a specific branch, you can specify this as well (here, `release` is used as
an example):
```
(@v1.8) pkg> add https://github.com/brainandforce/Xtal.jl#release
```
```julia-repl
julia> Pkg.add(url="https://github.com/brainandforce/Xtal.jl", rev="release")
```
The current development state is kept in the `main` branch, and the most recent stable version is
the head of the `release` branch. Specific releases for a minor version may be found by suffixing
the minor version with `/release`: for instance, `0.1/release` contains the latest release version
in the 0.1 series (currently 0.1.26).

## Current features

* Reading of common file formats:
     + abinit potential, density, and wavefunction outputs from versions 7.10.5 and 8.10.4
     + abinit HGH pseudopotentials
     + VASP POSCAR, WAVECAR, DOSCAR, and PROCAR
     + XCrysDen XSF
     + XYZ files
     + CPpackage2 outputs
* Writing of common file formats:
     + XCrysDen XSF
     + XYZ files
     + VASP POSCAR
     + LAMMPS atomic position data
* Operations on datagrids:
     + Addition, subtraction, multiplication
     + FFTs and inverse FFTs on real and reciprocal space data

## Planned features

This project is just starting to get off the ground, but here's what we have planned:

  * Reading and writing of common file formats:
      + XTL files
      + CIF files
      + abinit input and output files
  * Manipulation of data grids associated with crystal structures
      + Real space grid reinterpolation
  * Other operations
      + Converting k-point meshes to symmetry-reduced lists of weighted k-points
     
...and more that we will decide in time! If you'd like to contibute, be sure to read the included
[contributing guidelines.](CONTRIBUTING.md)

[docs-stable-img]:  https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]:  https://brainandforce.github.io/Xtal.jl/stable
[docs-dev-img]:     https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:     https://brainandforce.github.io/Xtal.jl/dev
[aqua-img]:         https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]:         https://github.com/JuliaTesting/Aqua.jl
