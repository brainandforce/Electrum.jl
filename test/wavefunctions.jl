wf = Electrum.readWAVECAR_new("files/WAVECAR")

@testset "Wavefunctions" begin
    i = PlanewaveIndex(1, 2, 3, 4, -5, 6)
    @test size(wf) === (1, 6, 12, 13, 13, 13)
    @test axes(wf) === (Base.OneTo(1), Base.OneTo(6), Base.OneTo(12), -6:6, -6:6, -6:6)
    wf[1, 2, 3, 4, -5, 6] = 420
    # Test that getindex and setindex! indices are the same
    @test wf[1,2,3,4,-5,6] == 420
    # Periodic indexing
    @test wf[1, 2, 3, 13, -13, 0] === wf[1, 2, 3, 0, 0, 0]
    @test_throws BoundsError wf[1, 7, 12, 0, 0, 0]
    @test wf[:,:,:] == wf[1, 1:6, 1:12] == wf
    @test wf[i] == wf[CartesianIndex(i)]
end
