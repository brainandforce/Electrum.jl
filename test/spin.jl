@testset "Spin" begin
    @test all(last(Multiplicity(n)) == (n - 1)//2 for n in 1:100)
    @test all(Multiplicity(n)[2] == (n - 1)//2 + 1 for n in 2:100)
    @test all(iszero(last(Multiplicity(n)) + first(Multiplicity(n))) for n in 1:100)
    @test all(Multiplicity(n) == UnitRange(Multiplicity(n)) for n in 1:100)
    @test all(
        UnitRange(Multiplicity(n)) === first(Multiplicity(n)):last(Multiplicity(n)) 
        for n in 1:100
    )
    @test all(Multiplicity(n)[2] == UnitRange(Multiplicity(n))[2] for n in 2:100)
    @test_throws AssertionError Multiplicity{0}()
    @test_throws AssertionError Multiplicity(0)
end