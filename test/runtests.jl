using LinearAlgebra, StaticArrays
using Test, Aqua, Electrum

Aqua.test_all(Electrum; project_toml_formatting=false)

xsf = readXSF3D("files/test.xsf")
v80_den = read_abinit_density("files/Sc_eq_o_DEN")
v80_wfk = read_abinit_wavefunction("files/Sc_eq_o_WFK")
# These two have the same data. Well...almost...
poscar = readPOSCAR("files/POSCAR")
lammps = read_lammps_data("files/lammps.data", atom_types = [14, 77])

@testset "All tests" begin
    include("internals.jl")
    include("lattices.jl")
    include("datagrids.jl")
    include("supercell.jl")
    include("crystals.jl")
    include("filetypes.jl")
    include("kpoints.jl")
    include("wavefunctions.jl")
end
