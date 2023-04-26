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

@testset "k-point mesh" begin
    kptmesh = KPointMesh(
        [
            KPoint(0.0, 0.0, 0.0, weight = 1),
            KPoint(0.2, 0.0, 0.0, weight = 6),
            KPoint(0.4, 0.0, 0.0, weight = 6),
            KPoint(0.2, 0.2, 0.0, weight = 6),
            KPoint(0.4, 0.2, 0.0, weight = 6),
            KPoint(0.0, 0.0, 1/3, weight = 2),
            KPoint(0.2, 0.0, 1/3, weight = 12),
            KPoint(0.4, 0.0, 1/3, weight = 12),
            KPoint(0.2, 0.2, 1/3, weight = 12),
            KPoint(0.4, 0.2, 1/3, weight = 12)
        ],
        diagm(SVector{3}(5,5,3))
    )
    # @test nkpt(kptmesh) == 10
    @test size(kptmesh) == (10,)
    @test axes(kptmesh) == (Base.OneTo(10),)
    @test sum(weight.(kptmesh)) == det(kptmesh.grid) == 75
    @test kptmesh[end] == KPoint(0.4, 0.2, 1/3, weight = 12)
    @test allunique(kptmesh)
    @test kptmesh[1:10] == kptmesh
end
