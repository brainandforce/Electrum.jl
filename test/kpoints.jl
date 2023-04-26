@testset "k-points" begin
    @test KPoint(0, 0, 0, weight = 1) != KPoint(0, 0, 0, weight = 2)
    @test hash(KPoint(0, 0, 0, weight = 1)) != hash(KPoint(0, 0, 0, weight = 2))
    @test KPoint(1, 2, 3) == KPoint(0, 0, 0)
    @test length(KPoint(1, 2, 3, 4, 5)) == 5
    @test length(KPoint{2}) == 2
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
    kptmesh[6] = [0.6, 0.9, 1/3]
    @test kptmesh[6] == KPoint{3}([0.6, 0.9, 1/3], 1)
    @test convert(KPointMesh, kptmesh.points) == KPointMesh(kptmesh.points)
    @test convert(KPointMesh, kptmesh.points) isa KPointMesh{3}
end
