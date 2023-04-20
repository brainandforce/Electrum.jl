# Fast Fourier transforms and related methods
# for `RealSpaceDataGrid`, `HKLData`, and `ReciprocalWavefunction`

"""
    fftfreq(g::RealSpaceDataGrid{D}) -> Array{SVector{D,Float64},D}

Returns the discrete Fourier transform frequency bins for a `RealSpaceDataGrid`. The frequency units
are angular wavenumbers, matching the 2π factors that are introduced when transforming between a
`RealBasis` and a `ReciprocalBasis`.
"""
function FFTW.fftfreq(g::RealSpaceDataGrid{D}) where D
    return collect(
        Iterators.product(
            (fftfreq(size(g)[d], 2π .* size(g)[d] / lengths(basis(g))[d]) for d in 1:D)...
        )
    )
end

"""
    fft(g::RealSpaceDataGrid) -> HKLData

Performs a fast Fourier transform on the data in a `RealSpaceDataGrid` and returns an `HKLData`.
"""
FFTW.fft(g::RealSpaceDataGrid) = HKLData(ReciprocalBasis(basis(g)), fft(g.data) * voxelsize(g))

"""
    ifft(g::HKLData) -> RealSpaceDataGrid

Performs an inverse fast Fourier transform on an `HKLData` and returns a `RealSpaceDataGrid`.
"""
FFTW.ifft(g::HKLData) = RealSpaceDataGrid(RealBasis(basis(g)), ifft(g.data) / voxelsize(g))
