@testset "Multiplicity" begin
    @test all(size(Multiplicity(n)) == (n,) for n in 1:100)
    @test all(last(Multiplicity(n)) == (n - 1)//2 for n in 1:100)
    @test all(Multiplicity(n)[2] == -(n - 1)//2 + 1 for n in 2:100)
    @test all(iszero(last(Multiplicity(n)) + first(Multiplicity(n))) for n in 1:100)
    @test all(Multiplicity(n) == UnitRange(Multiplicity(n)) for n in 1:100)
    @test all(
        UnitRange(Multiplicity(n)) === first(Multiplicity(n)):last(Multiplicity(n)) 
        for n in 1:100
    )
    @test all(Multiplicity(n)[2] == UnitRange(Multiplicity(n))[2] for n in 2:100)
    @test_throws BoundsError Multiplicity{3}()[0]
    @test_throws BoundsError Multiplicity{3}()[4]
    # Constructor robustness
    @test Multiplicity(1) === Multiplicity{1}()
    @test Multiplicity(1.0) === Multiplicity{1}()
    @test_throws AssertionError Multiplicity{0}()
    @test_throws AssertionError Multiplicity(0)
    @test_throws InexactError Multiplicity{1/2}()
    @test_throws InexactError Multiplicity(1/2)
    # Text representation
    @test eval(Meta.parse(repr(Multiplicity(3)))) === Multiplicity(3)
end

@testset "Spin bivectors" begin
    z_matrix = @SMatrix [0 1 0; -1 0 0; 0 0 0]
    s_z = SpinBivector(z_matrix)
    @test s_z[1,2] == 1
    @test s_z[2,1] == -1
    @test !Electrum._is_skew_symmetric([1 2 3; 4 5 6; 7 8 9])
    @test_throws AssertionError SpinBivector{3}([1 2 3; 4 5 6; 7 8 9])
    @test SpinBivector(SVector{3}(1, 0, 0), [0, 1, 0]) === s_z
    @test SpinBivector([0, 1, 0], SVector{3}(1, 0, 0)) === -s_z
    @test SpinBivector{3}([0 1 0; -1 0 0; 0 0 0]) === s_z
    @test SpinBivector{3,Int}([0 1 0; -1 0 0; 0 0 0]) === s_z
end
