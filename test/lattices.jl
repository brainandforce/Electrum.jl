@testset "Lattices" begin
    # Get basis vectors from the XSF file
    b = basis(xsf["this_is_3Dgrid#1"])
    # Inversion should give us 2π along the diagonal
    @test ReciprocalBasis(b) == ReciprocalBasis{3}(diagm([2π,2π,2π]) * Electrum.ANG2BOHR)
    # Check that conversion between real and reciprocal bases is invertible
    @test b ≈ RealBasis(ReciprocalBasis(b))
    # Dimensionality should be inferred when static matrices are used
    @test RealBasis(SMatrix{2,2}(1, 0, 0, 1)) isa RealBasis{2}
    # Wrong dimensionality should throw an exception
    @test_throws DimensionMismatch RealBasis{3}([1 0; 0 1])
    @test_throws DimensionMismatch convert(ReciprocalBasis{2}, b)
    # Gram matrix test
    @test gram(b) === (M = convert(SMatrix, b); M' * M)
end
