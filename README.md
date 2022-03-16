# Xtal.jl

A Julia package for working with crystal structures and various associated file formats, with the
aim of making theory development easier for everyone.

## How to install

Xtal.jl is not in the Julia package registry, so you'll need to install it like this:

```julia-repl
(@v1.7) pkg> add https://github.com/brainandforce/Xtal.jl 
```

You can access package mode by typing `]` at the REPL.

## Current features

* Reading of common file formats:
     + XCrysDen XSF data files
     + VASP and ABINIT outputs

## Planned features

This project is just starting to get off the ground, but here's what we have planned:

 * Reading and writing of common file formats:
     + XCrysDen XSF data files
     + VASP and ABINIT outputs
     + XYZ and XTL files
     + CIF files
 * Manipulation of data grids associated with crystal structures
     + Common elementwise mathematical operations on data grids
     + Conversion of real space data to reciprocal space data (and vice versa)
     
...and more that we will decide in time!