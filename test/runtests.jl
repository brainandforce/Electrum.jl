using LinearAlgebra, StaticArrays
using Test, Aqua, Xtal


xsf = readXSF3D("files/test.xsf")
v80_den = read_abinit_density("files/Sc_eq_o_DEN")
v80_wfk = read_abinit_wavefunction("files/Sc_eq_o_WFK")

@testset "Basis vectors" begin
    # Get basis vectors from the XSF file
    b = basis(xsf["this_is_3Dgrid#1"])
    # Inversion should give us 2π along the diagonal
    @test ReciprocalBasis(b) == ReciprocalBasis{3}(diagm([2π,2π,2π]))
    # Check that conversion between real and reciprocal bases is invertible
    @test b ≈ RealBasis(ReciprocalBasis(b))
end

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

@testset "XSF files" begin
    # Check that the key name is correct
    @test haskey(xsf.data, "this_is_3Dgrid#1")
    # Check that the size of the RealSpaceDataGrid is correct
    @test size(xsf["this_is_3Dgrid#1"]) == (4, 4, 4)
end

@testset "abinit outputs" begin
    den = v80_den["density_total"]
    wfk = v80_wfk["wavefunction"]
    # Check that the correct FFT grid size is read
    @test size(den) == (24, 24, 36)
    # Check that there's 1 spin, 4 k-points, and 8 bands
    @test size(wfk) == (1, 4, 8)
end
