# Fast Fourier transforms and related methods
# for `RealSpaceDataGrid`, `HKLData`, and `ReciprocalWavefunction`

"""
    fftfreq(g::RealSpaceDataGrid{D,<:Any}) -> Array{NTuple{D,Float64},D}

Returns the discrete Fourier transform frequency bins for a `RealSpaceDataGrid`.
"""
function FFTW.fftfreq(g::RealSpaceDataGrid{D,<:Any}) where D
    return collect(
        Iterators.product(
            (fftfreq(size(g)[d], size(g)[d] / lengths(basis(g))[d]) for d in 1:D)...
        )
    )
end

"""
    fft(g::RealSpaceDataGrid) -> HKLData

Performs a fast Fourier transform on the data in a `RealSpaceDataGrid` and returns an `HKLData`.
"""
function FFTW.fft(g::RealSpaceDataGrid)
    return HKLData(ReciprocalBasis(basis(g)), fft(grid(g)) * voxelsize(g))
end

"""
    ifft(g::HKLData) -> RealSpaceDataGrid

Performs an inverse fast Fourier transform on an `HKLData` and returns a `RealSpaceDataGrid`.
"""
function FFTW.ifft(hkl::HKLData)
    # TODO: Get voxel size for the associated RealSpaceDataGrid
    return RealSpaceDataGrid(RealBasis(basis(hkl)), ifft(grid(g)) / voxelsize(g))
end
