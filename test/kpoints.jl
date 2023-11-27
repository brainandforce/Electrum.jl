@testset "Shift vectors" begin
    # Equality and hashing
    @test KPoint(0, 0, 0, weight = 1) == KPoint(0, 0, 0, weight = 1)
    @test KPoint(0, 0, 0, weight = 1) != KPoint(0, 0, 0, weight = 2)
    @test hash(KPoint(0, 0, 0, weight = 1)) == hash(KPoint(0, 0, 0, weight = 1))
    @test hash(KPoint(0, 0, 0, weight = 1)) != hash(KPoint(0, 0, 0, weight = 2))
    @test zero(ShiftVector{Electrum.ByRealSpace,3}) != zero(KPoint{3})
    @test zero(KPoint{3,Float32}) == zero(KPoint{3,Int})
    # Constructors
    @test zero(KPoint{3}) === zero(KPoint{3,Bool})
    @test zero(KPoint{3}) === KPoint(false, false, false; weight = true)
    @test zero(KPoint{3,Int}) === KPoint(0, 0, 0, weight = 1)
    @test zero(KPoint{3,Int}) === KPoint(SVector{3}(0, 0, 0), 1)
    @test zero(KPoint{3,Float64}) === KPoint{3}(SVector{3}(0, 0, 0), 1.0)
    @test zero(KPoint{3,Float32}) === KPoint{3,Float32}(SVector{3}(0, 0, 0), 1.0)
    @test zero(KPoint{3,Float64}) === KPoint{3}([0, 0, 0], 1.0)
    @test zero(KPoint{3,Float32}) === KPoint{3,Float64}([0, 0, 0], 1.0)
    @test_throws Exception KPoint(SMatrix{3,1}([0, 0, 0]))
    @test_throws Exception KPoint{3}(SMatrix{3,1}([0, 0, 0]))
    @test_throws Exception KPoint{3,Int}(SMatrix{3,1}([0, 0, 0]))
    # Traits
    @test Electrum.DataSpace(zero(KPoint{3})) === Electrum.ByReciprocalSpace{3}()
    @test Electrum.ByCoordinate(zero(KPoint{3})) === Electrum.ByFractionalCoordinate{3}()
    # Truncation
    # TODO: see note for trunc() in src/vectors.jl
    @test truncate(KPoint(1, 2, 3)) == KPoint(0, 0, 0)
    @test truncate(KPoint(0.5, -0.5, 1.5)) == KPoint(0.5, -0.5, -0.5)
    @test truncate(KPoint(0.5, -0.5 + eps(Float64), 1.5)) == KPoint(0.5, -0.5 + eps(Float64), -0.5)
    @test truncate(KPoint(0.5, -0.5 - eps(Float64), 1.5)) == KPoint(0.5,  0.5 - eps(Float64), -0.5)
    # Length measurement
    @test length(KPoint(1, 2, 3, 4, 5)) == 5
    @test length(KPoint{2}) == 2
    @test size(KPoint(1, 2, 3, 4, 5)) == (5,)
    @test size(KPoint{2}) == (2,)
    @test convert(Vector, KPoint(0.1, 0.2, 0.3)) == [0.1, 0.2, 0.3]
    @test convert(Vector, KPoint(0.1, 0.2, 0.3)) isa Vector{<:Real}
    @test convert(SVector, KPoint(0.1, 0.2, 0.3)) === SVector{3}(0.1, 0.2, 0.3)
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
