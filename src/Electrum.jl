module Electrum

using LinearAlgebra
using StaticArrays
using FFTW
using Printf
using ComputedFieldTypes
using NormalForms

const ELEMENTS = 
( 
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "As", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
)

# TODO: can we modify this to make it unnecessary?
const ELEMENT_LOOKUP = Dict{String, Int}([(ELEMENTS[n] => n) for n in eachindex(ELEMENTS)])

# Matrices used to reduce conventional cells to corresponding primitive cells
# These can be accessed by calling the centering letter as a symbol
# So if you want the face-centered reduction matrix, use REDUCTION_MATRIX_3D[:F]
const REDUCTION_MATRIX_3D = 
(
    # Transpose these matrices in your head, please
    A = SMatrix{3,3,Float64}( 1,    0,    0,
                              0,  1/2,  1/2,
                              0, -1/2,  1/2),
    B = SMatrix{3,3,Float64}( 1/2,  0,  1/2,
                                0,  1,    0, 
                             -1/2,  0,  1/2),
    C = SMatrix{3,3,Float64}( 1/2,  1/2,  0,
                             -1/2,  1/2,  0,
                                0,    0,  1),
    F = SMatrix{3,3,Float64}( 1/2,  1/2,    0,
                              1/2,    0,  1/2,
                                0,  1/2,  1/2),
    I = SMatrix{3,3,Float64}(-1/2,  1/2,  1/2,
                              1/2, -1/2,  1/2, 
                              1/2,  1/2, -1/2),
    P = SMatrix{3,3,Float64}(LinearAlgebra.I)
)

"""
    Electrum.TOL_DEF

Default tolerance for discrepancies in floating point values.
"""
const TOL_DEF = 1e-8

"""
    Electrum.BOHR2ANG

Converts lengths in bohr to angstrom
"""
const BOHR2ANG = 0.52917720859

#= Units of these two constants should be energy^-1 * length^-2
"""
    Electrum.CABINIT

The value of 2m_e/ħ^2 as used in abinit (hartree^-1 * angstrom^-2).
"""
const CABINIT = 7.142129652186264
=#

"""
    Electrum.CVASP

The value of 2m_e/ħ^2 as used in VASP (eV^-1 * angstrom^-2, but not exactly).
"""
const CVASP = 0.262465831
# the correct value is actually 0.26246842360751754

# Functionality here has been superseded by allequal() in Julia 1.8.0
"""
    Electrum._allsame(itr)

Returns `true` if all the elements of an iterator are identical.
"""
_allsame(itr) = all(x -> x == first(itr), itr)

"""
    Electrum._selfdot(v)

Computes the dot product of a vector with itself.
"""
_selfdot(v) = dot(v,v)

# Methods used in vector operations that go beyond what's available in LinearAlgebra
include("math.jl")
# Abstract types used in type tree
include("types.jl")
export AbstractBasis, AbstractCrystal, AbstractCrystalData, AbstractRealSpaceData,
       AbstractReciprocalSpaceData, AbstractHKL, AbstractKPoints, AbstractDensityOfStates
# Methods and structs for working with crystal lattices
include("lattices.jl")
export RealBasis, ReciprocalBasis
export triangularize, lengths, volume, angles_cos, angles_rad, angles_deg, isdiag, qr, maxHKLindex
# Methods and structs for working with atomic positions
include("atoms.jl")
export NamedAtom, AbstractAtomPosition, FractionalAtomPosition, CartesianAtomPosition,
       AbstractAtomList, AtomList, PeriodicAtomList
export name, atomic_number, isdummy, displacement, occupancy, deduplicate, supercell, atomtypes,
       atomcounts, natomtypes
# Methods and structs for working with crystal data
include("crystals.jl")
export Crystal, CrystalWithDatasets
export data, generators, set_transform!
# Methods and structs for working with different types of data associated with crystals
include("data/realspace.jl")
export RealSpaceDataGrid
export basis, shift, grid, volume, voxelsize, coord, nearest_index, integrate
include("data/reciprocalspace.jl")
export KPointGrid, KPointList, BandAtKPoint, BandStructure, HKLData, HKLDict,
       ReciprocalWavefunction
export nkpt, nband, bounds
include("data/dos.jl")
export DensityOfStates, ProjectedDensityOfStates, FatBands
export fermi, smear, energies, nelectrons
include("data/atomic.jl")
export SphericalHarmonic
include("data/fft.jl")
export fftfreq, fft, ifft
# Methods and structs for working with different file formats
include("filetypes.jl")
export readXYZ, writeXYZ, readXSF3D, readXSF, writeXSF, readCPcoeff, readCPgeo, readCPcell
include("software/abinit.jl")
export read_abinit_density, read_abinit_potential, read_abinit_wavefunction, read_abinit_DEN,
       read_abinit_POT, read_abinit_anaddb_out, read_abinit_WFK, read_abinit_anaddb_in,
       write_abinit_modes, read_abinit_anaddb_PHDOS
include("software/vasp.jl")
export readPOSCAR, writePOSCAR4, readWAVECAR, readDOSCAR, readPROCAR, get_fermi, readKPOINTS
include("software/lammps.jl")
export write_lammps_data
# Show methods for pretty printing this module's structs
include("show.jl")
# Precompilation directives
include("precompile.jl")

end # end of module
