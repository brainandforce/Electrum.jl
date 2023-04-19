@testset "K-points" begin
    klist = KPointList(
        SVector{3,Float64}[
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5]
        ]
    )
    @test nkpt(klist) == 4
    @test length(klist) == 4
    @test size(klist) == (4,)
    @test axes(klist) == (Base.OneTo(4),)
    @test sum(klist.weights) ≈ 1
    # This will have to change when the weights are stored as integers
    # Direct comparison causes the weights to change since they are normalized to 1
    klist2 = KPointList(klist.points[1:3], klist.weights[1:3])
    @test klist[1:3] == klist2
    @test sum(klist.weights) ≈ 1
end
