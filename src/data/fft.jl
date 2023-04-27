# Fast Fourier transforms and related methods
# for `RealSpaceDataGrid`, `HKLData`, and `ReciprocalWavefunction`
"""
    fftfreq(g::AbstractDataGrid{D}) -> Array{SVector{D,Float64},D}

Returns the Fourier transform frequency bins for an `AbstractDataGrid`.

For real space data, the frequency bins will be angular wavenumbers, matching the 2π factors that
are introduced when transforming between a `RealBasis` and a `ReciprocalBasis`. The convention used
by `FFTW.fftfreq()` is also used: frequency bins at or above the Nyquist frequency will be negative.

For reciprocal space data, the frequencies are binned with the assumption that the lattice vectors
are given in angular wavenumbers, and they represent real space coordinates. The Nyquist frequency
convention is *not* used, so all elements will have positive indices.
"""
function FFTW.fftfreq(g::AbstractDataGrid, ::ByRealSpace)
    return map(i -> SVector(2π .* Tuple(i) ./ size(g)), FFTBins(g))
end

function FFTW.fftfreq(g::AbstractDataGrid, ::ByReciprocalSpace)
    return map(i -> SVector(Tuple(i) .* size(g) ./ 2π), CartesianIndices(g))
end

FFTW.fftfreq(g::AbstractDataGrid) = fftfreq(g, DataSpace(g))

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
