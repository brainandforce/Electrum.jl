@testset "Internals" begin
    # Type promotion
    @test Electrum.promote_typeof((1, Float16(2), Float32(3))...) === Float32
    # This actually works with both types and instances mixed together...
    @test Electrum.promote_eltype(Array{Float64}, zeros(SVector{3,Int})) === Float64
    # Linear independence
    @test Electrum.is_linearly_independent([1, 1], [2, 2]) === false
    # Conversion of scalars/vectors to transformation matrices
    import Electrum.convert_to_transform
    @test convert_to_transform(2, Val{3}()) == diagm(fill(2, SVector{3}))
    @test convert_to_transform(2, 3) == diagm(fill(2, 3))
    @test convert_to_transform(2*LinearAlgebra.I, 3) == diagm(fill(2, 3))
    @test convert_to_transform([1,2,3]) == diagm([1,2,3])
    @test convert_to_transform(SVector{3}(1,2,3)) == diagm(SVector{3}(1,2,3))
    @test convert_to_transform([1,2,3], Val{3}()) == diagm(SVector{3}(1,2,3))
    # Arguments not convertible to integers should fail
    @test_throws InexactError convert_to_transform(1.5, 2)
    # Data space traits
    @test BySpace(xsf["this_is_3Dgrid#1"]) === ByRealSpace()
    @test BySpace(v80_wfk["wavefunction"]) === ByReciprocalSpace()
    # Offset axes
    @test Base.has_offset_axes(xsf["this_is_3Dgrid#1"])
    @test Base.has_offset_axes(v80_wfk["wavefunction"])
    @test Base.has_offset_axes(v80_wfk["wavefunction"][1,1,1])
    # SUnitVector
    @test Electrum.SUnitVector{3}(2) == SVector{3,Bool}(0, 1, 0)
    @test Electrum.SUnitVector{3}(2) + Electrum.SUnitVector{3}(1) == SVector{3,Bool}(1, 1, 0)
    @test 2 * Electrum.SUnitVector{3}(3) === SVector{3}(0, 0, 2)
    @test_throws ErrorException Electrum.SUnitVector{3}(0)
    @test_throws ErrorException Electrum.SUnitVector{3}(4)
    @test Electrum.SUnitVector{3}(1)[1] === true
    @test Electrum.SUnitVector{3}(1)[:] === SVector{3,Bool}(1, 0, 0)
    @test Electrum.SUnitVector{3}(1)[1:2] == [1, 0]
    @test_throws MethodError Electrum.SUnitVector{3,Char}(1)
    # BySpace traits
    @test BySpace(RealBasis{3}) === ByRealSpace()
    @test BySpace(basis(wavecar)) === ByReciprocalSpace()
    @test BySpace(wavecar) === ByReciprocalSpace()
    @test BySpace(typeof(wavecar)) === ByReciprocalSpace()
    @test_throws Exception BySpace(420)
    @test_throws Exception BySpace(Int)
    @test inv(ByRealSpace) === ByReciprocalSpace
    @test inv(ByReciprocalSpace) === ByRealSpace
    @test inv(ByRealSpace()) === ByReciprocalSpace()
    @test inv(ByReciprocalSpace()) === ByRealSpace()
    # ByCoordinate traits
    @test ByCoordinate(KPoint(0, 0, 0)) === ByFractionalCoordinate()
    @test ByCoordinate(KPoint{3,Int}) === ByFractionalCoordinate()
    @test_throws Exception ByCoordinate(420)
    @test_throws Exception ByCoordinate(Int)
    # Requirements for compatible traits
    (k1, k2, s) = (KPoint(0, 0, 0), KPoint(1/2, 1/2, 1/2), zero(ShiftVector{ByRealSpace,3}))
    @test isnothing(Electrum.require_same_space(k1, k2))
    @test isnothing(Electrum.require_dual_space(k1, s))
    @test isnothing(Electrum.require_same_coordinate(k1, k2))
    @test_throws Exception Electrum.require_same_space(k1, s)
    @test_throws Exception Electrum.require_dual_space(k1, k2)
end

@testset "ZeroTo" begin
    # Constructors
    z = Electrum.ZeroTo(5)
    zz = Electrum.ZeroTo(UInt8(5))
    @test Electrum.ZeroTo{UInt8}(5) === zz
    @test_throws ArgumentError Electrum.ZeroTo{UInt8}(-4)
    # Equality
    @test z == zz
    # Element types
    @test eltype(z) == Int
    @test eltype(zz) == UInt8
    # Equivalence to UnitRange types
    @test z == 0:5
    @test zz == UInt8(0):UInt8(5)
    @test all(z .=== 0:5)
    @test first(z) === 0
    @test first(zz) === zero(UInt8)
    @test last(z) === 5
    @test last(zz) === UInt8(5)
    # Indexing
    @test z[1] === zero(eltype(z))
    @test z[end] === 5
    @test size(z) == tuple(6)
    @test length(z) == 6
    @test eachindex(z) === Base.OneTo(6)
    @test axes(z) === tuple(Base.OneTo(6))
    @test_throws BoundsError z[0]
    @test_throws BoundsError z[7]
end
