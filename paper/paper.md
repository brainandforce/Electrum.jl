---
title: 'Electrum.jl: a kit for building interactive chemical theory tools'
tags:
  - Julia
  - chemistry
  - density functional theory
  - crystallography
authors:
  - name: Brandon S. Flores
    orcid: 0000-0003-0816-4491
    affiliation: 1
    corresponding: true
  - name: Amber Lim
    orcid: 0000-0001-9893-1740
    affiliation: 1
  - name: Patrick K. Cross
    orcid: 0000-0002-4064-4033
    affiliation: 1
  - name: Leah C. Garman
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Michael E. Davies
    affiliation: 2
affiliations:
 - name: Department of Chemistry, Univesity of Wisconsin-Madison, USA
   index: 1
 - name: School of Computer, Data & Information Sciences, Univesity of Wisconsin-Madison, USA
date: 10 May 2023
bibliography: paper.bib

---

# Summary

Software implementing density functional theory (DFT) in a planewave basis allows for the
calculation of a large number of material properties - besides simply recovering the ground-state
wavefunction for the system, packages may implement tools for calculating phonon modes, electron
localization, and elasticity. However, structural chemists have developed their own theory tools for
understanding the factors that favor specific atomic arrangements: some examples of these tools are
crystal orbital Hamilton population (COHP) [@Blöchl:1993], which resolves chemical bonding by
interaction energies and atom type, and DFT-Chemical Pressure [@Fredrickson:2012], which projects
the steric effects in crystals onto atom-centered spherical harmonic representations of pressure.

Many tools utilized in structural chemistry are written in C or FORTRAN, which allow for high
performance - critical for systems with large numbers of degrees of freedom - but the low-level
nature of those languages serve as a barrier to entry for those who do not have a background in 
programming. In particular, we consider new chemistry graduate students whose first introduction
to software development comes with coursework or research. The Julia programming language can
achieve similar speeds to C or FORTRAN implementations while providing a lower barrier to entry for
novice programmers. With a toolkit that provides structures for a number of commonly used data types
and their associated methods, chemical theorists can focus on writing the most critical portions of
their theory tools without having to reimplement those structures for every new tool.

# Statement of need

`Electrum.jl` provides a number of data types that can be used to store and manipulate data
associated with crystal structures. This data may be read from a variety of sources: at the moment,
we support outputs from abinit [@Gonze:2020; @Romero:2020; @Gonze:2002], the Vienna *ab initio*
Simulation Package (VASP) [], and LAMMPS [], as well as the XCrysDen XSF data format [].

The data types provided by `Electrum.jl` are designed to allow for the storage of crystal data in
arbitrary dimensions. This functionality may be useful in the study of structural phenomena modeled
in projective spaces, such as incommensurately modulated structures and quasicrystals. 

`Electrum.jl` was designed for graduate and undergraduate students who are familiar with solid-state
chemistry but are new to programming. The Julia REPL provides a convenient interactive interface for
working with the data, and all exported data types have pretty-print methods that provide a summary
of relevant information in a data structure without polluting the output.

To maintain performance in a variety of scenarios, we took special care to optimize the layouts of
commonly used data structures. Static vectors and matrices from the `StaticArrays.jl` package are
used wherever possible in data structure definitions to match the dimensionality of the structure.

# Examples

## Constructing a supercell using arbitrary integer linear transformations

In molecular dynamics, it is often preferable to construct a simulation box that has orthogonal
basis vectors. Some crystal structures are not described with an orthogonal basis, but may be
transformed to one. This task is complicated by the need to generate new atomic positions, which is
straightforward for diagonal matrices, but significantly more complicated for non-diagonal matrices.

To facilitate this, the `supercell()` function allows for the conversion of a unit cell stored in a
`PeriodicAtomList` to a supercell transformed by an integer matrix, integer vector (corresponding to
a diagonal matrix), or integer scalar. As an example, one may transform data from a VASP POSCAR file
containing the basis vectors and atomic coordinates for MgZn2, a hexagonal Laves phase, and
constructs an orthogonal supercell.

```julia
using Electrum

MgZn2 = readPOSCAR("POSCAR")  # containing the atomic coordinates
T = [1 1 0; -1 1 0; 0 0 1] * []
MgZn2_sc = supercell(MgZn2, T)
writePOSCAR4("POSCAR2", MgZn2_sc)
```

## Fourier transforming DFT electron density

One can calculate a diffraction pattern using the results of a DFT calculation by performing a
Fourier transform of the electron density. The `fft()` and `ifft()` functions from `FFTW.jl` are
overloaded for use with `RealSpaceDataGrid` and `HKLData`, allowing for transformations between real
and reciprocal space representations of the same data. These Fourier transforms factor in the basis
vectors associated with the grids. 

A visual representation of the calculated diffraction pattern of {insert phase name} is given in
Figure ?. The calculated values were obtained from an abinit calculation {insert calculation details
here}. The phases associated with each reflection are encoded in the color of each spot. This is
compared to a real diffraction pattern obtained for the same compound using Mo-Kα radiation.
Discrepancies between calculated and real diffraction patterns may be attributed to the DFT
calculation being run at zero temperature and the approximations made by choice of the
exchange-correlation functional.

One may verify that the inverse FFT of the FFT of the density grid returns the same values to within
floating-point error:

# Development roadmap

`Electrum.jl` will seek to integrate itself with `AtomsBase.jl` to allow for easier interoperability
with other Julia packages that are used to perform chemical computation, such as `DFTK.jl`. This
integration will allow other packages to use methods provided by `Electrum.jl` to read files from
formats that were previously not supported, particularly the FORTRAN binary outputs from abinit and
VASP.

Data visualization was a primary consideration for the authors, and while visualization is outside
of the scope of `Electrum.jl`, we aim to provide `Makie.jl` plot recipes in the future, allowing for
convenient visualization of data structures within Visual Studio Code or other editors.

Further testing and development will also be directed at the use of GPUs to massively parallelize
large computations. This will leverage the data structures provided by packages such as `CUDA.jl`,
or `oneAPI.jl`.

# Acknowledgements

BSF, AL, PKC, and LCG all acknowledge Daniel C. Fredrickson for our training in solid-state and
computational chemistry, as well as discussions that were critical to many of the design decisions
made in the package.

BSF acknowledges Yingbo Ma and Chris Elrod for assistance in spinning off the
[HermiteNormalForms.jl](https://github.com/YingboMa/HermiteNormalForms.jl) package into the 
[NormalForms.jl](https://github.com/brainandforce/NormalForms.jl) package, which allows for the
computation of Smith and Hermite normal forms needed to perform supercell transformations. This
package is a dependency of `Electrum.jl`.

BSF, AL, and PKC acknowledge NSF grants (insert grant numbers for the Fredrickson group and the
department clusters).

# References
