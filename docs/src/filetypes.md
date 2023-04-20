# File formats

Electrum.jl supports a variety of different file formats from abinit, VASP, and LAMMPS.

# ABINIT

Outputs from abinit 7.10.5 and abinit 8.10.3 are supported. The headers from these versions contain
the numbers `57` and `80` respectively. Fundamentally, these files are FORTRAN binary files written
in sequential access mode with 4-byte record markers.

The file reading functionality from this package has not been tested on outputs from other versions
of abinit, but may work for other releases in versions 7 and 8 of abinit. However, these functions
will fail if the header does not contain either of those numbers in the header.

## Densities and potentials

Density outputs are suffixed with `_DEN` when written by abinit. Depending on the calculation type,
there can be either 1, 2, or 4 components (depending on the value of `nsppol`). The first component
is always the total density.

There is also the kinetic energy density with suffix `_KDEN` which should use the same format as
electron density files.

Potential outputs can have several different suffixes depending on the component of the potential
chosen to be written:
  * `_POT`: Total potential.
  * `_VPSP`: Local components of the pseudopotentials used.
  * `_VHA`: Hartree potential.
  * `_VHXC`: Sum of the Hartree and exchange-correlation potentials.
  * `_VXC`: Exchange-correlation potential. 

Both density and potential files follow the same format. The `read_abinit_density()` and 
`read_abinit_potential()` functions can be used to load in density and potential files.

## Wavefunctions

Wavefunctions can be read in by `read_abinit_wavefunction()`. This assumes that the wavefunction is
stored by k-points (`istwfk` should be equal to `nkpt*1` in the input file).

# VASP

Unfortunately, VASP file formats generally do not map neatly onto the data structures provided by
Electrum.jl.

## Files generated by VASP calculations

### DOSCAR

The DOSCAR file contains data need to plot a density of states curve.

### POSCAR and CONTCAR

The POSCAR file contains the basis vectors and all of the atomic positions used to generate a
crystal structure. 

The CONTCAR file is where the atomic positions from a calculation are written. In the case of
geometry optimizations, the contents of this file will differ.

### PROCAR

The PROCAR file contains data needed to plot fat bands (band structures that contain information
about contributions from atomic orbitals).

### WAVECAR

The WAVECAR file contains the coefficients for the wavefunction's reciprocal space representation
at each k-point. 

## Functions

### Reading

### Writing