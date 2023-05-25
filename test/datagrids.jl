@testset "DataGrid" begin
    g = xsf["this_is_3Dgrid#1"]
    # Equality and hashing checks
    h = deepcopy(g)
    @test g == h
    @test hash(g) == hash(h)
    h[0, 0, 0] = 420
    @test g != h
    @test hash(g) != hash(h)
    # Test broadcasting
    reference = RealDataGrid(2 * g.data, basis(g), shift(g))
    @test g + g == reference
    @test 2 * g == reference
    @test g .+ g == reference
    @test 2 .* g == reference
    @test volume(g) == Electrum.ANG2BOHR^3
    # similar() and zeros() constructors
    @test size(similar(g)) == size(g)
    @test eltype(similar(g)) == eltype(g)
    s = similar(g, Int, (4, 20, 69))
    @test s isa RealDataGrid{3,Int}
    @test size(s) === (4, 20, 69)
    @test basis(s) === basis(g)
    @test shift(s) === shift(g)
    z = zeros(ReciprocalDataGrid{3,Int}, basis(g), 4, 20, 69)
    @test z.data == zeros(Int, 4, 20, 69)
    @test basis(z) == ReciprocalBasis(g)
    @test shift(z) isa KPoint{3}
    @test weight(shift(z)) == 1
end

@testset "Fourier transforms" begin
    data = [cispi((x - 1)/5 + (y - 1)/10 + (z - 1)/15) for x in 1:10, y in 1:20, z in 1:30]
    g = RealDataGrid(data, RealBasis{3}(diagm([1,2,3])))
    h = fft(g)
    # Check that the Fourier transform and its inverse undo each other
    @test g ≈ ifft(fft(g))
    @test g ≈ fft(ifft(g))
    @test h ≈ fft(ifft(h))
    @test h ≈ ifft(fft(h)) 
    # Test FFT indexing mechanics
    @test collect(FFTLength(4)) == [0, 1, -2, -1]
    @test collect(FFTLength(5)) == [0, 1, 2, -2, -1]
    @test sort(FFTLength(4)) == -2:1
    @test sort(FFTLength(5)) == -2:2
    @test collect(FFTBins(4)) == [CartesianIndex(x) for x in [0, 1, -2, -1]]
    @test collect(FFTBins(3,3)) == [CartesianIndex(x,y) for x in [0, 1, -1], y in [0, 1, -1]]
    @test FFTBins(CartesianIndices(FFTBins(3,3))) == FFTBins(3,3)
end
