using LinearAlgebra, StaticArrays
using Test, Aqua, Electrum

tmpdir = mktempdir()

excluded_methods = Function[Base.unsafe_convert]
# Base.getindex ambiguity was resolved in 1.9:
# https://github.com/JuliaLang/julia/pull/41807
# TODO: can we still try to test other methods for Base.getindex? It's pretty important...
if VERSION < v"1.9"
    push!(excluded_methods, Base.getindex)
end
# Base.Sort.defalg ambiguity should be removed in 1.10:
# https://github.com/JuliaLang/julia/pull/47383
if VERSION < v"1.10"
    push!(excluded_methods, Base.Sort.defalg)
end

Aqua.test_all(Electrum;
    ambiguities = (exclude = excluded_methods,)
)

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
