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

    include("supercell.jl")

    @testset "RealSpaceDataGrid" begin
        g = xsf["this_is_3Dgrid#1"]
        # Check the indexing
        @test g[0] === 0.0
        @test g[1] === 1.000
        @test g[0,0,0] === 0.0
        @test g[1,1,1] === 1.732
        @test g[1,2,3] === 3.742
    end

    @testset "Fourier transforms" begin
        g = xsf["this_is_3Dgrid#1"]
        hkl = fft(g)
        # Check that the Fourier transform and its inverse undo each other
        @test g ≈ ifft(fft(g))
    end

    include("filetypes.jl")

end
