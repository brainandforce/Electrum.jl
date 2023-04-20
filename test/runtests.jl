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
    @testset "Basis vectors" begin
        # Get basis vectors from the XSF file
        b = basis(xsf["this_is_3Dgrid#1"])
        # Inversion should give us 2π along the diagonal
        @test ReciprocalBasis(b) == ReciprocalBasis{3}(diagm([2π,2π,2π]))
        # Check that conversion between real and reciprocal bases is invertible
        @test b ≈ RealBasis(ReciprocalBasis(b))
        # Dimensionality should be inferred when static matrices are used
        @test RealBasis(SMatrix{2,2}(1, 0, 0, 1)) isa RealBasis{2}
        # Wrong dimensionality should throw an exception
        @test_throws DimensionMismatch RealBasis{3}([1 0; 0 1])
        @test_throws DimensionMismatch convert(ReciprocalBasis{2}, b)
    end
    include("datagrids.jl")
    include("supercell.jl")
    include("crystals.jl")
    include("filetypes.jl")
    include("kpoints.jl")
    @testset "Miscellaneous" begin
        # Data space traits
        @test Electrum.data_space(xsf["this_is_3Dgrid#1"]) === Electrum.ByRealSpace{3}()
        @test Electrum.data_space(v80_wfk["wavefunction"]) === Electrum.ByReciprocalSpace{3}()
    end
end
