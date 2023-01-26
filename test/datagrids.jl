@testset "RealSpaceDataGrid" begin
    g = xsf["this_is_3Dgrid#1"]
    # Check the indexing
    @test g[0] === 0.0
    @test g[1] === 1.000
    @test g[0,0,0] === 0.0
    @test g[1,1,1] === 1.732
    @test g[1,2,3] === 3.742
end
@testset "Fourier transforms" begin
    g = xsf["this_is_3Dgrid#1"]
    hkl = fft(g)
    # Check that the Fourier transform and its inverse undo each other
    @test g ≈ ifft(fft(g))
end
