using LinearAlgebra, StaticArrays
using Test, Aqua, Electrum

Aqua.test_all(Electrum; project_toml_formatting=false)

xsf = readXSF3D("files/test.xsf")
v80_den = read_abinit_density("files/Sc_eq_o_DEN")
v80_wfk = read_abinit_wavefunction("files/Sc_eq_o_WFK")
poscar = readPOSCAR("files/")

@testset "All tests" begin
    @testset "Basis vectors" begin
        # Get basis vectors from the XSF file
        b = basis(xsf["this_is_3Dgrid#1"])
        # Inversion should give us 2π along the diagonal
        @test ReciprocalBasis(b) == ReciprocalBasis{3}(diagm([2π,2π,2π]))
        # Check that conversion between real and reciprocal bases is invertible
        @test b ≈ RealBasis(ReciprocalBasis(b))
    end
    include("datagrids.jl")
    include("supercell.jl")
    include("filetypes.jl")
end
