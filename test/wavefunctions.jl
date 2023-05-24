@testset "Wavefunctions" begin
    i = PlanewaveIndex(1, 2, 3, 4, -5, 6)
    @test size(wavecar) === (1, 6, 12, 13, 13, 13)
    @test axes(wavecar) === (Base.OneTo(1), Base.OneTo(6), Base.OneTo(12), -6:6, -6:6, -6:6)
    wavecar[1, 2, 3, 4, -5, 6] = 420
    @test wavecar[1, 1, 1, 4, 0, 0] == -4.900393f-5 + 1.3592809f-5im
    @test iszero(wavecar[1, 1, 1, 5, 0, 0])
    @test iszero(wavecar[1, 1, 1, 5, 0, 0])
    @test wavecar[1, 1, 1, -4, 0, 0] == -4.900393f-5 + 1.3592809f-5im
    @test wavecar[1, 6, 12, -1, -1, -1] == -0.070487864f0 - 0.1318744f0im
    # Test that getindex and setindex! indices are the same
    @test wavecar[1,2,3,4,-5,6] == 420
    # Periodic indexing
    @test wavecar[1, 2, 3, 13, -13, 0] === wavecar[1, 2, 3, 0, 0, 0]
    @test_throws BoundsError wavecar[1, 7, 12, 0, 0, 0]
    @test wavecar[:,:,:] == wavecar[1, 1:6, 1:12] == wavecar
    @test wavecar[i] == wavecar[CartesianIndex(i)]
    @test wavecar[1, 2, 3] isa ReciprocalDataGrid{3}
    @test vec(wavecar[1, 2, 3].data) == vec(wavecar.data[:, 3, 2, 1])
end
