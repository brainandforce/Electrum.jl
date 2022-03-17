"""
    Xtal

A module to assist with performing calculations on crystal structures.
"""
module Xtal

using   StaticArrays
using   LinearAlgebra
using   Printf

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
const ELEMENT_LOOKUP = Dict{String, Int}([(ELEMENTS[n] => n) for n in 1:length(ELEMENTS)])

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
    TOL_DEF

Default tolerance for discrepancies in floating point values.
"""
const TOL_DEF = 1e-8

const BOHR2ANG = 0.529177

# the value of 2m/ħ^2 as used in VASP
const CVASP = 0.262465831

"""
    _allsame(itr)

Returns `true` if all the elements of an iterator are identical.
"""
_allsame(itr) = all(x -> x == first(itr), itr)

"""
    _selfdot(v)

Computes the dot product of a vector with itself.
"""
_selfdot(v) = dot(v,v)

include("vectors.jl")
include("types.jl")
include("lattices.jl")
include("atoms.jl")
include("crystals.jl")
include("data.jl")
include("filetypes.jl")
include("show.jl")

# Abstract types to export
export  AbstractLattice, AbstractCrystal, AbstractCrystalData, AbstractRealSpaceData, 
        AbstractReciprocalSpaceData, AbstractKPoints, AbstractDensityOfStates
# Concrete types to export
export  RealLattice, ReciprocalLattice, Crystal, CrystalWithDatasets, RealSpaceDataGrid, 
        KPointGrid, KPointList, HKLData, ReciprocalWavefunction, DensityOfStates,
        ProjectedDensityOfStates, BandAtKPoint, BandStructure, AtomList, AtomPosition
# Functions to export
export  nkpt, nband, lattice2D, lattice3D, readXSF3D, readWAVECAR, read_abinit_density,
        readDOSCAR, data, basis, gridsize, writeXSF, coord, atomname, atomicno

#= For debugging purposes only
export  read_abinit_header, read_abinit_header_57, read_abinit_header_80, get_abinit_version
=#

end # end of module
