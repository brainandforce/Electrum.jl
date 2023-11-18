@testset "Lattices" begin
    # Get basis vectors from the XSF file - the basis associated with the datagrid
    b = basis(xsf["this_is_3Dgrid#1"])
    # Inversion should give us 2π along the diagonal
    @test ReciprocalBasis(b) ≈ ReciprocalBasis{3}(diagm(fill(2π / Electrum.ANG2BOHR, 3)))
    # Check that conversion between real and reciprocal bases is invertible
    @test b ≈ RealBasis(ReciprocalBasis(b))
    # Dimensionality should be inferred when static matrices are used
    @test RealBasis(SMatrix{2,2}(1, 0, 0, 1)) isa RealBasis{2}
    # Wrong dimensionality should throw an exception
    @test_throws DimensionMismatch RealBasis{3}([1 0; 0 1])
    @test_throws DimensionMismatch convert(ReciprocalBasis{2}, b)
    # Gram matrix test
    @test gram(b) === (M = convert(SMatrix, b); M' * M)
    # Type promotion
    @test promote_type(RealBasis{3,Int64}, RealBasis{3,Float32}) === RealBasis{3,Float32}
    @test promote_type(RealBasis{3,Int64}, ReciprocalBasis{3,Float32}) === SMatrix{3,3,Float32,9}
    # Metrics
    @test all(x == Electrum.ANG2BOHR for x in lengths(b))
    # eachvertex iterator
    @test isequal(
        Electrum._gen_vertices(reshape(1:69, 23, 3), 1:420),
        Electrum._gen_vertices(1:69, reshape(1:420, 3, 35, 4))
    )
    @test collect(eachvertex([6 0; 0 9])) == vec(collect(eachvertex(SMatrix{2,2}(6, 0, 0, 1))))
    result = SVector{3,Float64}[
        [0, 0, 0], [1, 0, 0], [0, 2, 0], [1, 2, 0],
        [0, 0, 3], [1, 0, 3], [0, 2, 3], [1, 2, 3]
    ]
    @test all(isequal.(eachvertex(diagm([1, 2, 3])), result))
    @test isequal(
        collect(eachvertex(SMatrix{3,3}(diagm([1, 2, 3])))),
        reshape(result, 2, 2, 2)
    )
end
