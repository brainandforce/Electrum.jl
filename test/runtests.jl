using LinearAlgebra, StaticArrays
using Test, Xtal

xsf = readXSF3D("files/test.xsf")

@testset "Basis vectors" begin
    # Get basis vectors from the XSF file
    real_basis = basis(xsf.gen)
    # Check that conversion between real and reciprocal bases is invertible
    @test (RealBasis(ReciprocalBasis(b)).vs .- b.vs) ≈ 0
end

@testset "File formats" begin
    @test haskey(xsf.data, "this_is_3Dgrid#1")
    @test size(xsf["this_is_3Dgrid#1"]) == (4, 4, 4)
end
