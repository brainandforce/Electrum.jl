# Xtal.jl

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

If you'd like to track a specific branch, you can specify this as well (here, `dev` is used as an
example):

```
(@v1.8) pkg> add https://github.com/brainandforce/Xtal.jl#dev
```
```julia-repl
julia> Pkg.add(url="https://github.com/brainandforce/Xtal.jl", rev="dev")
```

Our current branching strategy involves the `main` branch (tracked by default), the `dev` branch 
(contains features and fixes that will make it to the next patch to the current minor version), and
the `next` branch (contains potentially breaking changes to the next minor version).

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
     + FFTs on real space data grids

## Planned features

This project is just starting to get off the ground, but here's what we have planned:

 * Reading and writing of common file formats:
     + XTL files
     + CIF files
     + abinit input and output files
 * Manipulation of data grids associated with crystal structures
     + Real space grid reinterpolation
     + Gradient calculations
     
...and more that we will decide in time! If you'd like to contibute, be sure to read the included
[contributing guidelines.](CONTRIBUTING.md)
