module Electrum

using LinearAlgebra
using StaticArrays
using FFTW
using Printf
using InlineStrings
using ComputedFieldTypes
using NormalForms
using Requires

import InlineStrings.InlineString15

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
    Electrum.BOHR2ANG

Converts lengths in bohr to angstrom.
"""
const BOHR2ANG = 0.52917721090380

"""
    Electrum.ANG2BOHR

Converts lengths in angstrom to bohr.
"""
const ANG2BOHR = 1.8897261246229133

"""
    Electrum.HARTREE2EV

Converts energies in hartree to electron-volts.
"""
const HARTREE2EV = 27.21138624598853

"""
    Electrum.EV2HARTREE

Converts energies in electron-volts to hartree.
"""
const EV2HARTREE = 0.03674932217565428

"""
    Electrum.CVASP

The value of 2m_e/ħ^2 as used in VASP (eV^-1 * angstrom^-2, but not exactly).
"""
const CVASP = 0.262465831
# the correct value is actually 0.26246842360751754

#---Convenience methods----------------------------------------------------------------------------#

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

"""
    Electrum.hash_fields(x, h::UInt)

Hashes all fields of a data structure recursively.
"""
hash_fields(x, h::UInt) = (for n in fieldnames(typeof(x)); h = hash(getfield(x, n), h); end; h)

#---Type promotion---------------------------------------------------------------------------------#
"""
    Electrum.promote_typeof(args...)

Calls `promote_type()` on the types of the arguments, equivalent to `promote_type(typeof.(x)...)`.

This function is implemented in Julia Base, but it is not part of the public API; a stable
implementation is provided here.
"""
promote_typeof(x...) = promote_type(typeof.(x)...)

"""
    Electrum.promote_eltype(args...)

Calls `promote_type()` on the element types of the arguments, equivalent to 
`promote_type(eltype.(x)...)`.

This function is implemented in Julia Base, but it is not part of the public API; a stable
implementation is provided here.
"""
promote_eltype(x...) = promote_type(eltype.(x)...)

#---Included files---------------------------------------------------------------------------------#

# Methods used in array operations that go beyond what's available in LinearAlgebra
include("math.jl")
export FFTBins, FFTLength
# Dispatch traits for data
include("traits.jl")
export BySpace, ByRealSpace, ByReciprocalSpace, ByCoordinate, ByCartesianCoordinate, 
    ByFractionalCoordinate
# Coordinate vectors
include("geometry/vectors.jl")
export CoordinateVector, RealCartesianCoordinate, RealFractionalCoordinate, 
    ReciprocalCartesianCoordinate, ReciprocalFractionalCoordinate, ShiftVector, KPoint
export weight
# Methods and structs for working with crystal lattices
include("geometry/lattices.jl")
export LatticeBasis, RealBasis, ReciprocalBasis
export eachvertex, basis, dual, dualbasis, lengths, volume, angles_cos, angles_rad, angles_deg, 
    gram, isdiag, qr, triangularize, maxHKLindex
# Map array elements to unit cell geometry
include("geometry/arraymaps.jl")
export LatticeDataMap
# Methods and structs for working with atomic positions
include("atoms.jl")
export NamedAtom, AbstractAtomPosition, FractionalAtomPosition, CartesianAtomPosition,
    AbstractAtomList, AtomList, PeriodicAtomList
export name, atomic_number, isdummy, displacement, occupancy, distance, deduplicate, supercell,
    atomtypes, atomcounts, natomtypes
# Methods and structs for working with crystal data
include("crystals.jl")
export AbstractCrystal, Crystal, CrystalWithDatasets
export data, generators, get_transform, set_transform!
# Weighed k-points and k-point meshes
include("data/kpoints.jl")
export KPointMesh
export nkpt
# Energy/occupancy pairs
include("data/energies.jl")
export AbstractEnergyData, EnergyOccupancy, StateDensity, EnergiesOccupancies
export energy, occupancy, density, energies, occupancies, densities, min_energy, max_energy,
    min_occupancy, max_occupancy, fermi
# Spin data: multiplicity and direction
include("data/spin.jl")
export Multiplicity, SpinBivector
# Real and reciprocal space data grids
include("data/grids.jl")
export DataGrid, RealDataGrid, ReciprocalDataGrid
export shift, fft, ifft, fftfreq, voxelsize, integrate, partial_derivative, cell_gradient, gradient,
    directional_derivative, remove_shift
# New data grid implementation
include("data/newgrids.jl")
export LatticeDataMap
# Planewave wavefunctions
include("data/wavefunctions.jl")
export PlanewaveIndex, PlanewaveWavefunction
export nspin, nband, nonzero_g_indices, nonzero_g_vectors
# Band structures
include("data/bands.jl")
export BandAtKPoint, BandStructure
export nband
# Density of states
include("data/dos.jl")
export AbstractDensityOfStates, DensityOfStates, ProjectedDensityOfStates, FatBands
export smear, energies, nelectrons
# Data associated with atoms
include("data/atomic.jl")
export SphericalHarmonic
# Methods and structs for working with different file formats
include("filetypes.jl")
export readXYZ, writeXYZ
include("software/xcrysden.jl")
export readXSF3D, readXSF, writeXSF
include("software/cppackage.jl")
export readCPcoeff, readCPgeo, readCPcell
include("software/abinit.jl")
export read_abinit_DEN, read_abinit_POT, read_abinit_WFK
#= Removed exports for abinit anaddb functionality
export read_abinit_anaddb_out, read_abinit_anaddb_in, write_abinit_modes, read_abinit_anaddb_PHDOS
=#
include("software/vasp.jl")
export readPOSCAR, readCONTCAR, writePOSCAR, writeCONTCAR, readWAVECAR, readDOSCAR, readPROCAR,
    get_fermi
include("software/lammps.jl")
export read_lammps_data, write_lammps_data
# Show methods for pretty printing this module's structs
include("show.jl")
# Precompilation directives
include("precompile.jl")

function __init__()
    @require TOML="fa267f1f-6049-4f14-aa54-33bafae1ed76" include("software/toml.jl")
end

end # end of module
