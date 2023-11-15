@testset "Spin" begin
    @test all(last(SpinRange(n)) == (n - 1)//2 for n in 1:100)
    @test all(SpinRange(n)[2] == (n - 1)//2 + 1 for n in 2:100)
    @test all(iszero(last(SpinRange(n)) + first(SpinRange(n))) for n in 1:100)
    @test all(SpinRange(n) == UnitRange(SpinRange(n)) for n in 1:100)
    @test all(UnitRange(SpinRange(n)) === first(SpinRange(n)):last(SpinRange(n)) for n in 1:100)
    @test all(SpinRange(n)[2] == UnitRange(SpinRange(n))[2] for n in 2:100)
    @test_throws AssertionError SpinRange{0}()
    @test_throws AssertionError SpinRange(0)
end
