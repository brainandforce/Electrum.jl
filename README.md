# Xtal.jl

A Julia package for working with crystal structures, associated data, and various file formats,
with the aim of making theory development easier for everyone.

## How to install

Xtal.jl is not in the Julia package registry yet, so you'll need to install it like this:

```julia-repl
(@v1.7) pkg> add https://github.com/brainandforce/Xtal.jl 
```

You can access package mode by typing `]` at the REPL.

## Current features

* Reading of common file formats:
     + abinit potential, density, and wavefunction outputs from versions 7.10.5 and 8.10.4
     + abinit HGH pseudopotentials
     + VASP WAVECAR, DOSCAR, and PROCAR
     + XCrysDen XSF
     + XYZ files
     + CPpackage2 outputs
* Writing of common file formats:
     + XCrysDen XSF
     + XYZ files
* Operations on datagrids:
     + Addition, subtraction, multiplication
     + FFTs on real space data grids

## Planned features

This project is just starting to get off the ground, but here's what we have planned:

 * Reading and writing of common file formats:
     + XTL files
     + CIF files
 * Manipulation of data grids associated with crystal structures
     + Real space grid reinterpolation
     
...and more that we will decide in time! If you'd like to contibute, be sure to read the included
contributing guidelines.