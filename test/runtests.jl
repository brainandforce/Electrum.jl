using LinearAlgebra, StaticArrays
using Test, Aqua, Electrum

tmpdir = mktempdir()

Aqua.test_all(Electrum; project_toml_formatting=false)

xsf = readXSF3D("files/test.xsf")
header_den = Electrum.read_abinit_header("files/Sc_eq_o_DEN")
header_wfk = Electrum.read_abinit_header("files/Sc_eq_o_WFK")
v80_den = read_abinit_DEN("files/Sc_eq_o_DEN")
v80_wfk = read_abinit_WFK("files/Sc_eq_o_WFK")
wavecar = readWAVECAR("files/WAVECAR")
# These two have the same data. Well...almost...
poscar = readPOSCAR("files/POSCAR")
lammps = read_lammps_data("files/lammps.data", atom_types = [14, 77])

@testset "All tests" begin
    include("internals.jl")
    include("lattices.jl")
    include("atoms.jl")
    include("datagrids.jl")
    include("supercell.jl")
    include("crystals.jl")
    include("filetypes.jl")
    include("kpoints.jl")
    include("wavefunctions.jl")
end
