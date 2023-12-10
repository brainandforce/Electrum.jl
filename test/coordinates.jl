@testset "Coordinates" begin
    c1 = RealCartesianCoordinate(1, 2, 3)
    # Constructors
    @test RealCartesianCoordinate(true, UInt8(2), 3) === c1
    @test RealCartesianCoordinate((true, UInt8(2), 3)) === c1
    @test RealCartesianCoordinate(MVector(1, 2, 3)) === c1
    @test_throws Exception RealCartesianCoordinate(SMatrix{1,3}(1, 2, 3))
    @test RealCartesianCoordinate{3}((1, 2, 3)) === c1
    @test RealCartesianCoordinate{3}([1, 2, 3]) === c1
    @test RealCartesianCoordinate{3}(MVector(1, 2, 3)) === c1
    @test RealCartesianCoordinate{3}(SMatrix{1,3}(1, 2, 3)) === c1
    @test RealCartesianCoordinate{3,Int}((1, 2, 3)) === c1
    @test RealCartesianCoordinate{3,Int}([1, 2, 3]) === c1
    @test RealCartesianCoordinate{3,Int}(MVector(1, 2, 3)) === c1
    @test RealCartesianCoordinate{3,Int}(SMatrix{1,3}(1, 2, 3)) === c1
    @test CoordinateVector{ByRealSpace}(c1) === c1
    @test CoordinateVector(c1) === c1
    # Conversion
    @test convert(SVector, c1) === SVector(1, 2, 3)
    @test convert(SVector{3}, c1) === SVector(1, 2, 3)
    @test convert(SVector{3,Float32}, c1) === SVector{3,Float32}(1, 2, 3)
    @test convert(MVector, c1) == MVector(1, 2, 3)
    c1_f32 = RealCartesianCoordinate{3,Float32}(1, 2, 3)
    @test convert(RealCartesianCoordinate{3,Float32}, c1) === c1_f32
    @test convert.(Float32, c1) === c1_f32
    @test_throws Exception convert(RealFractionalCoordinate{3}, c1)
    @test_throws Exception convert(ReciprocalCartesianCoordinate{3}, c1)
    # Zero vectors
    @test zero(RealFractionalCoordinate{3,Int}) === RealFractionalCoordinate(0, 0, 0)
    @test zero(RealFractionalCoordinate{3}) == RealFractionalCoordinate(0, 0, 0)
    # Math operations should maintain the correct type
    @test -c1 === RealCartesianCoordinate(-1, -2, -3)
    @test c1 + RealCartesianCoordinate(4, 5, 6) === RealCartesianCoordinate(5, 7, 9)
    @test RealCartesianCoordinate(5, 7, 9) - RealCartesianCoordinate(4, 5, 6) === c1
    @test c1 + SVector(4, 5, 6) === RealCartesianCoordinate(5, 7, 9)
    @test c1 + [4, 5, 6] === RealCartesianCoordinate(5, 7, 9)
    @test_throws Exception c1 + RealFractionalCoordinate(4, 5, 6)
    @test_throws Exception c1 + ReciprocalCartesianCoordinate(4, 5, 6)
    # Arbitrary multiplication operations do not retain information about space/coordinate system.
    @test c1 * 2 === SVector(2, 4, 6)
    @test 2 * c1 === SVector(2, 4, 6)
    @test c1 / 2 === SVector(1/2, 2/2, 3/2)
    @test 2 / c1 === SVector(2/1, 2/2, 2/3)
    @test 2 \ c1 === SVector(1/2, 2/2, 3/2)
    @test c1 \ 2 === SVector(2/1, 2/2, 2/3)
    @test dot(c1, c1) === dot(SVector(c1), SVector(c1)) === 14
    # Show methods should produce valid Julia code
    @test eval(Meta.parse(repr(c1))) == c1
    @test eval(Meta.parse(repr(RealFractionalCoordinate(1.0, 2.0, 3.0)))) == convert.(Float64, c1)
end

@testset "Shift vectors" begin
    # Equality and hashing
    @test KPoint(0, 0, 0, weight = 1) == KPoint(0, 0, 0, weight = 1)
    @test KPoint(0, 0, 0, weight = 1) != KPoint(0, 0, 0, weight = 2)
    @test hash(KPoint(0, 0, 0, weight = 1)) == hash(KPoint(0, 0, 0, weight = 1))
    @test hash(KPoint(0, 0, 0, weight = 1)) != hash(KPoint(0, 0, 0, weight = 2))
    @test zero(ShiftVector{Electrum.ByRealSpace,3}) != zero(KPoint{3})
    @test zero(KPoint{3,Float32}) == zero(KPoint{3,Int})
    # Constructors and zero k-point
    @test zero(KPoint{3}) === zero(KPoint{3,Bool})
    @test zero(KPoint{3}) === KPoint(false, false, false; weight = true)
    @test zero(KPoint{3,Int}) === KPoint(0, 0, 0, weight = 1)
    @test zero(KPoint{3,Int}) === KPoint(SVector{3}(0, 0, 0), 1)
    @test zero(KPoint{3,Float64}) === KPoint{3}(SVector{3}(0, 0, 0), 1.0)
    @test zero(KPoint{3,Float32}) === KPoint{3,Float32}(SVector{3}(0, 0, 0), 1.0)
    @test zero(KPoint{3,Float64}) === KPoint{3}([0, 0, 0], 1.0)
    @test zero(KPoint{3,Float32}) === KPoint{3,Float32}([0, 0, 0], 1.0)
    @test_throws Exception KPoint(SMatrix{3,1}([0, 0, 0]))
    @test_throws InexactError KPoint{3,Int8}(420, 69, 1337)
    # Traits
    @test BySpace(zero(KPoint{3})) === ByReciprocalSpace()
    @test ByCoordinate(zero(KPoint{3})) === ByFractionalCoordinate()
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
    @test convert(ReciprocalFractionalCoordinate, KPoint(0.1, 0.2, 0.3)) ===
        ReciprocalFractionalCoordinate(0.1, 0.2, 0.3)
    @test convert(KPoint, ReciprocalFractionalCoordinate(1, 2, 3)) === KPoint(1, 2, 3)
    @test convert(KPoint{3,Float64}, ReciprocalFractionalCoordinate(1, 2, 3)) === 
        KPoint(1.0, 2.0, 3.0)
    @test_throws Exception convert(KPoint, ReciprocalCartesianCoordinate(1, 2, 3))
    @test_throws Exception convert(KPoint, RealFractionalCoordinate(1, 2, 3))
    # Show methods should produce valid Julia code
    @test eval(Meta.parse(repr(KPoint(0.1, 0.2, 0.3)))) === KPoint(0.1, 0.2, 0.3)
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
