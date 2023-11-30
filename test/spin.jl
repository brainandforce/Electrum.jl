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
    x = [1, 0, 0]
    y = [0, 1, 0]
    (sx, sy) = SVector{3}.((x, y))
    z_matrix = @SMatrix [0 1 0; -1 0 0; 0 0 0]
    s_z = SpinBivector(z_matrix)
    @test Electrum._is_skew_symmetric(z_matrix)
    @test !Electrum._is_skew_symmetric([1 2 3; 4 5 6; 7 8 9])
    @test_throws AssertionError SpinBivector{3}([1 2 3; 4 5 6; 7 8 9])
    @test SpinBivector(sx, y) === s_z
    @test SpinBivector(y, sx) === -s_z
    @test SpinBivector{3}([0 1 0; -1 0 0; 0 0 0]) === s_z
    @test SpinBivector{3,Int}([0 1 0; -1 0 0; 0 0 0]) === s_z
    # Constructor with wedge product
    @test SpinBivector(sx, sy) === s_z
    @test SpinBivector{3}(x, y) === s_z
    @test SpinBivector{3,Float64}(x, y) == s_z
    # Wedge product kernel for AbstractVector
    @test Electrum._wedge_matrix(sx, sy) == Electrum._wedge_matrix(x, y)
    # Properties
    @test :data in propertynames(s_z, private = true)
    @test only(propertynames(s_z)) == :matrix
    # Indexing
    @test s_z[2,1] == -1
    @test s_z[1,2] == 1
    @test s_z[2] == -1
    @test s_z[4] == 1
    # Bad constructor inputs
    @test_throws ErrorException SpinBivector(x, y)
    @test_throws ErrorException SpinBivector(@SVector [1, 2, 3])
    @test_throws ErrorException SpinBivector{3}(@SVector [1, 2, 3])
    @test_throws ErrorException SpinBivector{3,Int}(@SVector [1, 2, 3])
    @test_throws ErrorException SpinBivector(SArray{Tuple{2,2,2}}(1, 2, 3, 4, 5, 6, 7, 8))
    @test_throws DimensionMismatch SpinBivector(SVector{2}(1, 0), SVector{3}(0, 1, 0))
    @test_throws DimensionMismatch SpinBivector{3}([1, 0, 0], [0, 1])
    # Traits
    @test Electrum.DataSpace(SpinBivector{3}) === Electrum.ByRealSpace{3}()
    @test Electrum.DataSpace(s_z) === Electrum.ByRealSpace{3}()
end
