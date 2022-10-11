using LinearAlgebra, StaticArrays
using Test, Xtal

xsf = readXSF3D("files/test.xsf")

@testset "Basis vectors" begin
    # Get basis vectors from the XSF file
    b = basis(xsf)
    # Check that conversion between real and reciprocal bases is invertible
    # The difference between their elements should be very small (though not exactly zero)
    @test all(RealBasis(ReciprocalBasis(b)).vs[m][n] ≈ b.vs[m][n] for m in 1:3, n in 1:3)
end

@testset "Fourier transforms" begin
    g = xsf["this_is_3Dgrid#1"]
    @test 
end

@testset "File formats" begin
    # Check that the key name is correct
    @test haskey(xsf.data, "this_is_3Dgrid#1")
    # Check that the size of the RealSpaceDataGrid is correct
    @test size(xsf["this_is_3Dgrid#1"]) == (4, 4, 4)
end
