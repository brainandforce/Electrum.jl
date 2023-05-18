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
   index: 2
date: 18 May 2023
bibliography: paper.bib

---

# Summary

Software implementing density functional theory (DFT) in a planewave basis allows for the
calculation of a large number of material properties - besides simply recovering the ground-state
wavefunction for the system, packages may implement tools for calculating phonon modes, electron
localization, and elasticity. However, structural chemists have developed their own tools to
understand the factors that favor specific atomic arrangements: some examples are crystal orbital
Hamilton population (COHP) [@Dronskowski:1993], which resolves chemical bonding by interaction
energies and atom type, and DFT-Chemical Pressure [@Fredrickson:2012], which projects the steric
effects in crystals onto atom-centered spherical harmonic representations of pressure.

Many tools utilized in structural chemistry are written in C or FORTRAN, which allow for high
performance - critical for systems with large numbers of degrees of freedom - but the low-level
nature of those languages serve as a barrier to entry for those who do not have a background in 
programming. In particular, we consider new chemistry graduate students who are introduced to
software development through coursework or research. The Julia programming language [@Bezanson:2017] 
can achieve similar speeds to C or FORTRAN implementations while providing a lower barrier to entry
for novice programmers. With a toolkit that provides structures for a number of commonly used data
types and their associated methods, chemical theorists can write extensible theory tools which can
easily interface with other Julia packages or computational chemistry software.

# Statement of need

`Electrum.jl` was designed for graduate and undergraduate students who may be introduced to
programming through computational chemistry. The package provides a number of data types that can be
used to store data associated with crystal structures, and methods that perform common operations.
As much as possible, we accomplish this by overloading Julia `Base` methods so that syntax is
intuitive for new users. Users may import data from a variety of sources: at the moment, we support
FORTRAN binary outputs from abinit [@Gonze:2020; @Romero:2020; @Gonze:2002] and Vienna *ab initio*
Simulation Package (VASP) [@Kresse:1993; @Kresse:1996A; @Kresse:1996B], atomic position data from
LAMMPS [@Thompson:2022], and the XCrysDen XSF data format [@Kokalj:1999]. Atomic position data may
be exported to the XYZ and VASP POSCAR formats, and more complex scalar data may be exported to the
XCrysDen XSF format.

Although `Electrum.jl` aims for simplicity for new users, the data types provide a great deal of 
flexibility for the construction of advanced tools. When possible, data structures are defined with
a dimension parameter that may be set arbitrarily, and methods which operate on them are defined
generically so that they produce valid results regardless of dimension. This enables the development
of tools that model crystals with more complex periodicity, such as incommensurately modulated
crystals and quasicrystals, whose lattices are commonly modeled in a higher-dimensional superspace 
which is cut at an irrational angle and projected down to 3 dimensions [@deBruijn:1981].

# Examples

## Constructing a supercell using arbitrary integer linear transformations

In molecular dynamics, it is often preferable to construct a simulation box that has orthogonal
basis vectors. Some crystal structures are not described with an orthogonal basis, but may be
transformed to one. This task is complicated by the need to generate new atomic positions, which is
straightforward for diagonal matrices, but significantly more complicated for non-diagonal matrices.

To facilitate this, the `supercell()` function allows for the conversion of a unit cell stored in a
`PeriodicAtomList` to a supercell transformed by an integer matrix, integer vector (corresponding to
a diagonal matrix), or integer scalar. As an example, one may transform data from a VASP POSCAR file
containing the basis vectors and atomic coordinates for MgZn₂, a hexagonal Laves phase, and
constructs an orthogonal supercell.

```julia
using Electrum

MgZn2 = readPOSCAR("POSCAR")  # containing the atomic coordinates
T = [1 1 0; -1 1 0; 0 0 1] * diagm([5, 5, 3])
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
calculation being run at zero temperature and the approximations made by choice of pseudopotential
and exchange-correlation functional.

One may verify that the inverse FFT of the FFT of the density grid returns the same values to within
floating-point error:

# Development roadmap

`Electrum.jl` seeks to integrate itself with `AtomsBase.jl` to allow for easier interoperability
with other Julia packages used to perform chemical computation, such as `DFTK.jl`. This integration
will allow other packages to use methods provided by `Electrum.jl` to read files in formats not 
supported by those packages, particularly the FORTRAN binary outputs from abinit and VASP.

Data visualization was a primary consideration for the authors, and while visualization is outside
of the scope of `Electrum.jl`, we aim to provide `Makie.jl` plot recipes in the future as part of a
supplemental package, allowing for convenient visualization of data structures in editors that can
render HTML.

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
