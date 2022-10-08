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
    fft(g::RealSpaceDataGrid{D,<:Number}; maxhkl=zeros(Int,D)) -> HKLData{D,<:Complex}

Performs a fast Fourier transform on the data in a `RealSpaceDataGrid{D,<:Number}` and
generates an `HKLData{D,<:Complex}`.
"""
function FFTW.fft(
    g::RealSpaceDataGrid{D,<:Number};
    maxhkl::AbstractVector{<:Integer} = zeros(Int, D)
) where D 
    # Generate the bounds needed for the HKLdata
    bounds = [range(div(sz, 2, RoundUp) - sz, length = sz) for sz in gridsize(g)]
    @debug string("Bounds:\t", bounds)
    # Calculate the shift factors to put the values in the right places
    shifts = [last(b) for b in bounds]
    # Take the grid fft
    # Permute the elements to get the indexing working
    # Multiply by the size of a voxel to get the right scaling
    f = circshift(fft(grid(g)), shifts .+ 1) .* voxelsize(g)
    # With defined bounds, truncate the data
    if all(!iszero, maxhkl)
        ind = [(-x:x) .- first(b) .+ 1 for (x,b) in zip(abs.(maxhkl), bounds)]
        f = f[ind...]
        bounds = [-x:x for x in abs.(maxhkl)]
    end
    return HKLData(f, bounds)
end
