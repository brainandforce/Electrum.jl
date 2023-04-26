@testset "RealSpaceDataGrid" begin
    g = xsf["this_is_3Dgrid#1"]
    # Check the indexing
    @test g[0] === 0.0
    @test g[1] === 1.000
    @test g[0,0,0] === 0.0
    @test g[1,1,1] === 1.732
    @test g[1,2,3] === 3.742
    # Equality and hashing checks
    h = deepcopy(g)
    @test g == h
    @test hash(g) == hash(h)
    h[0, 0, 0] = 420
    @test g != h
    @test hash(g) != hash(h)
end

@testset "Fourier transforms" begin
    g = xsf["this_is_3Dgrid#1"]
    hkl = fft(g)
    # Check that the Fourier transform and its inverse undo each other
    @test g ≈ ifft(fft(g))
    # Test FFT indexing mechanics
    @test collect(FFTLength(4)) == [0, 1, -2, -1]
    @test collect(FFTLength(5)) == [0, 1, 2, -2, -1]
    @test sort(FFTLength(4)) == -2:1
    @test sort(FFTLength(5)) == -2:2
    @test collect(FFTBins(4)) == [CartesianIndex(x) for x in [0, 1, -2, -1]]
    @test collect(FFTBins(3,3)) == [CartesianIndex(x,y) for x in [0, 1, -1], y in [0, 1, -1]]
    @test FFTBins(CartesianIndices(FFTBins(3,3))) == FFTBins(3,3)
end
