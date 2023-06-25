@testset "Internals" begin
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
    @test Electrum.DataSpace(xsf["this_is_3Dgrid#1"]) === Electrum.ByRealSpace{3}()
    @test Electrum.DataSpace(v80_wfk["wavefunction"]) === Electrum.ByReciprocalSpace{3}()
    # Offset axes
    @test Base.has_offset_axes(xsf["this_is_3Dgrid#1"])
    @test Base.has_offset_axes(v80_wfk["wavefunction"])
    @test Base.has_offset_axes(v80_wfk["wavefunction"][1,1,1])
    # SUnitVector
    @test Electrum.SUnitVector{3}(2) === SVector{3,Bool}(0, 1, 0)
    @test Electrum.SUnitVector{3}(2) + Electrum.SUnitVector{3}(1) == SVector{3,Bool}(1, 1, 0)
    @test 2 * Electrum.SUnitVector{3}(3) === SVector{3}(0, 0, 2)
    @test_throws ErrorException Electrum.SUnitVector{3}(0)
    @test_throws ErrorException Electrum.SUnitVector{3}(4)
    @test Electrum.SUnitVector{3}(1)[1] === true
    @test Electrum.SUnitVector{3}(1)[:] === SVector{3,Bool}(1, 0, 0)
    @test Electrum.SUnitVector{3}(1)[1:2] == [1, 0]
end
