@testset "Supercell construction" begin
    # Coordinates for Laves-MgZn2
    MgZn2 = PeriodicAtomList(
        RealBasis(5.223 * SMatrix{3,3,Float64}(1, 0, 0, -1/2, sqrt(3/4), 0, 0, 0, sqrt(8/3))),
        [
            FractionalAtomPosition(12, SVector(1/3, 2/3, 1/16)),
            FractionalAtomPosition(12, SVector(1/3, 2/3, 7/16)),
            FractionalAtomPosition(12, SVector(2/3, 1/3, 9/16)),
            FractionalAtomPosition(12, SVector(2/3, 1/3, 15/16)),
            FractionalAtomPosition(30, SVector(0.0, 0.0, 0.0)),
            FractionalAtomPosition(30, SVector(0.0, 0.0, 1/2)),
            FractionalAtomPosition(30, SVector(5/6, 2/3, 1/4)),
            FractionalAtomPosition(30, SVector(5/6, 1/6, 1/4)),
            FractionalAtomPosition(30, SVector(1/3, 1/6, 1/4)),
            FractionalAtomPosition(30, SVector(1/6, 1/3, 3/4)),
            FractionalAtomPosition(30, SVector(1/6, 5/6, 3/4)),
            FractionalAtomPosition(30, SVector(2/3, 5/6, 3/4)),
            # Two duplicates for good measure!
            FractionalAtomPosition(12, SVector(4/3, 2/3, 1/16)),
            FractionalAtomPosition(30, SVector(0.0, 0.0, 1/2))
        ]
    )
    # The extra Zn atom should be stripped because it's identical to an earlier one
    @test length(MgZn2) == 13
    # Build an orthorhombic supercell (treat hexagonal as orthorhombic C)
    transform = [1 1 0; -1 1 0; 0 0 1]
    MgZn2_sc = supercell(MgZn2, transform)
    # Check that the basis is actually orthorhombic
    @test convert(SMatrix, basis(MgZn2_sc)) ≈ diagm(SVector(9.0465013679, 5.223, 8.5291232843))
    # Check that the atomic list size scales with the transform matrix
    @test length(MgZn2_sc)/length(MgZn2) <= abs(det(transform))
    # Check that the duplicate Mg atom was stripped
    @test length(MgZn2_sc) == 24
    # Check each individual atom:
    #=
    @testset "MgZn2 supercell construction check" begin
        # Just work with the coordinates
        Mg = [coord(a) for a in filter(12, MgZn2_sc)]
        Zn = [coord(a) for a in filter(30, MgZn2_sc)]
        @test any(isapprox(SVector(1/2,   1/6,  1/16)), Mg)
        @test any(isapprox(SVector(1/2,   1/6,  7/16)), Mg)
        @test any(isapprox(SVector(1/2,   5/6,  9/16)), Mg)
        @test any(isapprox(SVector(1/2,   5/6, 15/16)), Mg)
        @test any(isapprox(SVector(  0,   2/3,  1/16)), Mg)
        @test any(isapprox(SVector(  0,   2/3,  7/16)), Mg)
        @test any(isapprox(SVector(  0,   1/3,  9/16)), Mg)
        @test any(isapprox(SVector(  0,   1/3, 15/16)), Mg)
        @test any(isapprox(SVector(0.0,   0.0,   0.0)), Zn)
        @test any(isapprox(SVector(  0,     0,   1/2)), Zn)
        @test any(isapprox(SVector(1/2,   1/2,   0.0)), Zn)
        @test any(isapprox(SVector(1/2,   1/2,   1/2)), Zn)
        @test any(isapprox(SVector(  0,   1/6,   1/4)), Zn)
        @test any(isapprox(SVector(  0,   5/6,   3/4)), Zn)
        @test any(isapprox(SVector(1/4,  5/12,   1/4)), Zn)
        @test any(isapprox(SVector(3/4,  5/12,   1/4)), Zn)
        @test any(isapprox(SVector(1/2,   1/3,   1/4)), Zn)
        @test any(isapprox(SVector(1/4, 11/12,   1/4)), Zn)
        @test any(isapprox(SVector(3/4, 11/12,   1/4)), Zn)
        @test any(isapprox(SVector(1/4,  1/12,   3/4)), Zn)
        @test any(isapprox(SVector(3/4,  1/12,   3/4)), Zn)
        @test any(isapprox(SVector(1/2,   1/3,   3/4)), Zn)
        @test any(isapprox(SVector(1/4,  7/12,   3/4)), Zn)
        @test any(isapprox(SVector(3/4,  7/12,   3/4)), Zn)
    end
    =#
end
