"""
    Electrum.PlanewaveIndex{D}

A special indexing type used to index the components of wavefunctions in a planewave basis.

In many computational chemistry packages, the standard indexing of wavefunction components occurs in
the following canonical order: spins, k-points, bands, then the h, k, and l indices of the G-vectors
associated with the coefficients. However, in many cases, users will want to select a spin, k-point,
or band before selecting a G-vector index.

To keep the syntax intuitive while maintaining performance (mostly by ensuring that all of the
components of a wavefunction are stored compactly) this index type ensures that the iteration occurs
in a natural order for users and in an efficient order for Julia.
    
This also ensures that G-vectors with negative indices are handled correctly: while the canonical
G-vectors are those within some defined ranges of indices, out of bounds indices are reinterpreted
automatically using modulo arithmetic.
"""
struct PlanewaveIndex{D}
    spin::Int
    kpoint::Int
    band::Int
    g::CartesianIndex{D}
    PlanewaveIndex(spin, kpt, band, g::CartesianIndex{D}) where D = new{D}(spin, kpt, band, g)
end

PlanewaveIndex(spin, kpt, band, g...) = PlanewaveIndex(spin, kpt, band, CartesianIndex(g))
PlanewaveIndex(spin, kpt, band, g::NTuple) = PlanewaveIndex(spin, kpt, band, CartesianIndex(g))
PlanewaveIndex(spin, kpt, band, g::StaticVector) = PlanewaveIndex(spin, kpt, band, Tuple(g))
PlanewaveIndex(i::CartesianIndex) = PlanewaveIndex(Tuple(i)...)

function PlanewaveIndex{D}(spin, kpt, band, gvec::AbstractVector) where D
    return PlanewaveIndex(spin, kpt, band, ntuple(i -> gvec[i], Val{D}()))
end

Base.show(io::IO, i::PlanewaveIndex) = print(io, PlanewaveIndex, (i.spin, i.kpoint, i.band, i.g))

Tuple(i::PlanewaveIndex) = (i.spin, i.kpoint, i.band, Tuple(i.g)...)
CartesianIndex(i::PlanewaveIndex) = CartesianIndex(Tuple(i))

Base.convert(T::Type{<:Tuple}, i::PlanewaveIndex) = Tuple(i)::T
Base.convert(T::Type{<:CartesianIndex}, i::PlanewaveIndex) = CartesianIndex(i)::T
Base.convert(T::Type{<:PlanewaveIndex}, i::CartesianIndex) = PlanewaveIndex(i)::T

#=
"""
    Electrum.PlanewaveIndices{D} <: AbstractArray{PlanewaveIndex{D},D}

An array of valid indices of a `PlanewaveWavefunction{D,T}` within the G-vector bounds determined by
the energy cutoff of the calculation that generated the data.
"""
struct PlanewaveIndices{D} <: AbstractArray{PlanewaveIndex{D},D}
    spins::Base.OneTo{Int}
    kpoints::Base.OneTo{Int}
    bands::Base.OneTo{Int}
    grange::NTuple{D,UnitRange{Int}}
    function PlanewaveIndices(
        spins::AbstractUnitRange,
        kpoints::AbstractUnitRange,
        bands::AbstractUnitRange,
        grange::NTuple{D,<:AbstractUnitRange}
    ) where D
        return new{D}(spins, kpoints, bands, grange)
    end
end

function PlanewaveIndices(
    spins::AbstractUnitRange,
    kpts::AbstractUnitRange,
    bands::AbstractUnitRange,
    gs::AbstractUnitRange...
)
    return PlanewaveIndices(spins, kpts, bands, gs)
end

function PlanewaveIndices(spins::Integer, kpts::Integer, bands::Integer, gs::AbstractUnitRange...)
    return PlanewaveIndices(Base.OneTo(spins), Base.OneTo(kpts), Base.OneTo(bands), gs)
end

Base.axes(p::PlanewaveIndices) = (p.spins, p.kpoints, p.bands, p.grange...)
Base.size(p::PlanewaveIndices) = length.(axes(p))

function Base.show(io::IO, p::PlanewaveIndices)
    print(
        io, PlanewaveIndices,
        (last(p.spins), last(p.kpoints), last(p.bands), p.grange...)
    )
end

Base.show(io::IO, ::MIME"text/plain", p::PlanewaveIndices) = show(io, p)

# Convert to a vector within the allowed range of G-vectors
# In other words, favor a negative index in that range over a larger positive index
function Base.getindex(p::PlanewaveIndices{D}, i::PlanewaveIndex{D}) where D
    i.spin in p.spins || throw(BoundsError(p, i))
    i.kpoint in p.kpoints || throw(BoundsError(p, i))
    i.band in p.bands || throw(BoundsError(p, i))
    sz = length.(p.grange)
    c = CartesianIndex(mod.(Tuple(i.g) .+ div.(sz, 2), sz) .- div.(sz, 2))
    return PlanewaveIndex(i.spin, i.kpoint, i.band, c)
end

Base.getindex(p::PlanewaveIndices, ::Integer, ::Integer, ::Integer) = CartesianIndices(p.grange)

function Base.getindex(
    p::PlanewaveIndices{D},
    spin::Integer,
    kpt::Integer,
    band::Integer,
    i::Vararg{<:Integer,D}    
) where D
    return p[PlanewaveIndex(spin, kpt, band, i...)]
end

Base.LinearIndices(p::PlanewaveIndices) = LinearIndices((p.grange..., p.bands, p.kpoints, p.spins))
=#

"""
    PlanewaveWavefunction{D,T} <: AbstractDataGrid{D,T}

Stores the components of a wavefunction constructed from a planewave basis. Usually, the coefficient
data type `T` will be a `ComplexF32`, as in DFT calculations, double precision convergence of the
density will correspond to single-precision converegnce of the wavefunction.

Internally, coefficients are stored in an `Array{4,T}`. Indexing is then manually implemented, with
a `D`-dimensional `CartesianIndex` used for accessing each coefficient associated with a G-vector.
`PlanewaveWavefunction` instances are mutable, with `getindex()` and `setindex!()` defined for them,
but they are not resizable, and the backing array should not be resized.
"""
struct PlanewaveWavefunction{D,T} <: AbstractDataGrid{D,T}
    basis::ReciprocalBasis{D}
    spins::Vector{SVector{D,Float64}}
    kpoints::KPointMesh{D}
    energies::Array{Float64,3}
    occupancies::Array{Float64,3}
    grange::NTuple{D,UnitRange{Int}}
    data::Array{T,4}
    function PlanewaveWavefunction(
        basis::AbstractBasis{D},
        spins::AbstractVector{<:StaticVector{D,<:Real}},
        kpoints::AbstractVector{KPoint{D}},
        energies::AbstractArray{<:Real,3},
        occupancies::AbstractArray{<:Real,3},
        grange::NTuple{D,<:AbstractUnitRange{<:Integer}},
        data::Array{T,4}
    ) where {D,T}
        @assert length(spins) === size(data, 4) "Mismatch in the number of spins"
        @assert length(kpoints) === size(data, 3) "Mismatch in the number of k-points"
        @assert size(energies) === size(data)[2:4] "Mismatch in the size of energy data"
        @assert size(occupancies) === size(data)[2:4] "Mismatch in the size of occupancy data"
        @assert prod(length.(grange)) === size(data, 1) "G-vector limits do not match data size"
        return new{D,T}(basis, spins, kpoints, energies, occupancies, grange, data)
    end
end

"""
    PlanewaveWavefunction{D,T}(
        basis::AbstractBasis{D},
        nspin::Integer,
        nkpt::Integer,
        nband::Integer,
        grange::AbstractUnitRange{<:Integer}...
    )

Constructs an empty `PlanewaveWavefunction` with `nspin` spins, `nkpt` k-points, `nband` bands, and
G-vectors in the ranges given by `grange`.
"""
function PlanewaveWavefunction{D,T}(
    basis::AbstractBasis,
    nspin::Integer,
    nkpt::Integer,
    nband::Integer,
    grange::Vararg{AbstractUnitRange{<:Integer},D}
) where {D,T}
    return PlanewaveWavefunction(
        basis,
        zeros(SVector{D,Float64}, nspin),
        KPointMesh(zeros(KPoint{D}, nkpt)),
        zeros(Float64, nband, nkpt, nspin),
        zeros(Float64, nband, nkpt, nspin),
        grange,
        zeros(T, prod(length.(grange)), nband, nkpt, nspin)
    )
end

Base.size(wf::PlanewaveWavefunction) = (reverse(size(wf.energies))..., length.(wf.grange)...)
Base.axes(wf::PlanewaveWavefunction) = (reverse(axes(wf.energies))..., wf.grange...)

# PlanewaveIndices(wf::PlanewaveWavefunction) = PlanewaveIndices(axes(wf)...)

function Base.LinearIndices(wf::PlanewaveWavefunction)
    return LinearIndices((wf.grange..., wf.bands, wf.kpoints, wf.spins))
end

nspin(wf::PlanewaveWavefunction) = length(wf.spins)
nkpt(wf::PlanewaveWavefunction) = length(wf.kpoints)
nband(wf::PlanewaveWavefunction) = size(wf.energies, 1)

# Override some of the more generic AbstractDataGrid methods
Base.getindex(wf::PlanewaveWavefunction, i...) = throw(MethodError(getindex, (wf, i...)))
Base.setindex!(wf::PlanewaveWavefunction, i...) = throw(MethodError(setindex!, (wf, i...)))

function Base.getindex(wf::PlanewaveWavefunction{D}, i::PlanewaveIndex{D}) where D
    all(in.((i.spin, i.kpoint, i.band), axes(wf)[1:3])) || throw(BoundsError(wf, Tuple(i)))
    l = LinearIndices(wf.grange)[CartesianIndex(mod.(Tuple(i.g), length.(wf.grange)) .+ 1)]
    @inbounds return wf.data[l, i.band, i.kpoint, i.spin]
end

Base.getindex(wf::PlanewaveWavefunction, i::CartesianIndex) = wf[PlanewaveIndex(i)]

function Base.getindex(wf::PlanewaveWavefunction, spin, kpt, band)
    # Convert integers to ranges so that arrays are sliced into smaller arrays
    s = spin isa Integer ? (spin:spin) : spin
    k = kpt isa Integer ? (kpt:kpt) : kpt
    b = band isa Integer ? (band:band) : band
    return PlanewaveWavefunction(
        wf.basis,
        wf.spins[s],
        wf.kpoints[k],
        wf.energies[b, k, s],
        wf.occupancies[b, k, s],
        wf.grange,
        wf.data[:, b, k, s]
    )
end

function Base.getindex(wf::PlanewaveWavefunction, spin::Integer, kpt::Integer, band::Integer)
    return HKLData(
        basis(wf),
        reshape(wf.data[:, band, kpt, spin], length.(wf.grange)),
        wf.kpoints[kpt]
    )
end

function Base.getindex(
    wf::PlanewaveWavefunction{D},
    spin::Integer,
    kpt::Integer,
    band::Integer,
    g::Vararg{<:Integer,D}
) where D
    return wf[PlanewaveIndex(spin, kpt, band, g)]
end

function Base.setindex!(wf::PlanewaveWavefunction{D}, x, i::PlanewaveIndex{D}) where D
    l = LinearIndices(wf.grange)[CartesianIndex(mod.(Tuple(i.g), length.(wf.grange)) .+ 1)]
    wf.data[l, i.band, i.kpoint, i.spin] = x
end

function Base.setindex!(
    wf::PlanewaveWavefunction{D},
    x,
    spin::Integer,
    kpt::Integer,
    band::Integer,
    g::Vararg{<:Integer,D}
) where D
    wf[PlanewaveIndex(spin, kpt, band, g)] = x
end

#---New WAVECAR and WFK reading functions----------------------------------------------------------#
"""
    readWAVECAR_new(file; quiet = false) -> PlanewaveWavefunction{3,Float32}

Reads a WAVECAR file output from a VASP 4.6 calcuation to the new `PlanewaveWavefunction` type.

Information about VASP WAVECAR files and much of the code was pulled from the WaveTrans website
(originally written in FORTRAN): https://www.andrew.cmu.edu/user/feenstra/wavetrans/

This function is limited to WAVECAR files which have an RTAG value of 45200 (meaning the data is
given as a `Complex{Float32}`) and have only a collinear magnetic field applied, like WaveTrans. It
should also be noted that the weights of the k-points are not present in the WAVECAR file, and are
set to 1 by default.

By default, the function is verbose, with output printed for every k-point parsed, due to the large
size of the wavefunction files. If this behavior is undesirable, the `quiet` keyword argument may be
set to `true`.
"""
function readWAVECAR_new(io::IO; quiet = false)
    # Function to increment HKL values in place 
    function incrementHKL!(hkl::AbstractVector{<:Integer}, bounds::AbstractVector{<:AbstractRange})
        # Loop through the vector indices, but in most cases we don't need them all
        for n in eachindex(hkl)
            # Increment the current vector component
            hkl[n] = (hkl[n] + 1 in bounds[n] ? hkl[n] + 1 : minimum(bounds[n]))
            # Only increment the next components if the current one is zero
            hkl[n] == 0 || break
        end
    end
    # Data entry counter (for the entries in the WAVECAR)
    count = 0
    # Number of bytes per record
    nrecl = Int(read(io, Float64))
    @debug "Record length: " * string(nrecl)
    # Number of spin components
    nspin = Int(read(io, Float64))
    # Check for the proper format
    rtag = Int(read(io, Float64))
    rtag == 45200 || error("Unsupported format: format value is " * string(rtag))
    # Jump to the next record
    count +=1; seek(io, count*nrecl)
    # Number of k-points
    nkpt = Int(read(io, Float64))
    # Number of bands
    nband = Int(read(io, Float64))
    # Energy cutoff
    ecut = read(io, Float64)
    # Reciprocal lattice vectors
    rlatt = convert(ReciprocalBasis, RealBasis{3}([read(io, Float64) for a = 1:3, b = 1:3]))
    # Get HKL coefficient bounds (as done in WaveTrans)
    hklbounds = SVector{3,UnitRange{Int}}(-g:g for g in maxHKLindex(rlatt, ecut))
    # Bare wavefunction to be filled
    wf = PlanewaveWavefunction{3,Complex{Float32}}(rlatt, nspin, nkpt, nband, hklbounds...)
    # Loop through the spins
    for s in 1:nspin
        # Loop through the k-points
        for kp in 1:nkpt
            # Seek to the next data
            count += 1; seek(io, count*nrecl)
            # Number of plane waves for this k-point
            pos = position(io)
            @debug string("File pointer at ", pos, " (", count, " * ", nrecl, ")")
            npw = Int(read(io, Float64))
            # Add the position of the k-point to the list
            wf.kpoints[kp] = [read(io, Float64) for n in 1:3]
            # Get energies and occupancies
            for b in 1:nband
                # Ordering is reversed (see PlanewaveIndex above...)
                wf.energies[b, kp, s] = read(io, Float64)
                skip(io, 8)
                wf.occupancies[b, kp, s] = read(io, Float64)
            end
            quiet || @info string(
                "Read in data for k-point ", kp, "/", nkpt, " (", npw, " planewaves/band)\n",
                "Reciprocal space coordinates: ", @sprintf("[%f %f %f]", wf.kpoints[kp]...)
            )
            for b in 1:nband
                # Seek to the next entry
                count +=1; seek(io, count*nrecl)
                # Reset the HKL indices
                hkl = zeros(MVector{3,Int})
                for p in 1:npw
                    # Get the planewave component
                    pw = read(io, Complex{Float32})
                    # Increment the HKL indices
                    while true
                        # Get the energy of the vector
                        sumkg = [dot(wf.kpoints[kp] + hkl, rlatt[:,n]) for n in 1:3]
                        etot = _selfdot(sumkg)/CVASP
                        # Break when the G-vector energy is below ecut
                        # This may occur immediately if the k-vector already meets the criteria
                        etot < ecut ? break : incrementHKL!(hkl, hklbounds)
                    end
                    # Store the data at the HKL index
                    # Note: data is stored first by k-points, then by bands
                    wf[s, kp, b, hkl...] = pw
                    # Increment it for the next iteration
                    incrementHKL!(hkl, hklbounds)
                end
            end
        end
    end
    return wf
end

readWAVECAR_new(filename; quiet = false) = open(io -> readWAVECAR_new(io; quiet), filename)

function read_abinit_wavefunction_new(io::IO)
    # Get the header from the file
    header = read_abinit_header(io)
    # Get the reciprocal lattice
    rlatt = convert(ReciprocalBasis, RealBasis(header.rprimd))
    # Get the minimum and maximum HKL values needed
    # Units for c (2*m_e/ħ^2) are hartree^-1 bohr^-2
    # this should affect only the size of the preallocated arrays
    bounds = SVector{3,UnitRange{Int}}(
        -g:g for g in maxHKLindex(rlatt, header.ecut, c=2)
    )
    @debug string(
        "hklbounds: ", bounds, "\n",
        "Calculated from ecut = ", header.ecut, " Hartree"
    )
    nb = maximum(header.nband)
    wf = PlanewaveWavefunction{3,Complex{Float64}}(rlatt, header.nsppol, header.nkpt, nb, bounds...)
    # Loop over all the spin polarizations
    for sppol in 1:header.nsppol
        # Loop over all the k-points
        for kpt in 1:header.nkpt
            # Skip over record length
            read(io, Int32)
            # Write the number of plane waves, spinors, and bands
            (npw, nspinor, nband) = [read(io, Int32) for n in 1:3]
            # Skip over another pair of record lengths
            read(io, Int32); read(io, Int32)
            # HKL indices of the band components
            hklinds = [read(io, SVector{3,Int32}) for n in 1:npw]
            # Skip over another pair of record lengths
            read(io, Int32); read(io, Int32)
            # Eigenvalues and band occupancy at the current k-point
            wf.energies[:, kpt, sppol] = [read(io, Float64) for n in (1:nband)]
            wf.occupancies[:, kpt, sppol] = [read(io, Float64) for n in (1:nband)]
            read(io, Int32)
            # Loop over all the bands (given in previous entry, not from header)
            for band in 1:nband
                read(io, Int32)
                cg = [read(io, Complex{Float64}) for m in 1:npw, n in 1:nspinor]
                # Now iterate through the supplied indices and set them
                for (ind, coeff) in zip(hklinds, cg)
                    wf[sppol, kpt, band, ind...] = coeff
                end
                # Skip marker
                read(io, Int32)
            end
            @info string(
                "Read in data for k-point $kpt/", header.nkpt, " ($npw planewaves/band)\n",
                "Reciprocal space coordinates: ", @sprintf("[%f %f %f]", header.kpt[kpt]...)
            )
        end
    end
    return CrystalWithDatasets{3,String,PlanewaveWavefunction{3,Complex{Float64}}}(
        Crystal(header), 
        Dict{String,PlanewaveWavefunction}("wavefunction" => wf)
    )
end

"""
    read_abinit_wavefunction_new(file)
        -> CrystalWithDatasets{3,String,PlanewaveWavefunction{3,Complex{Float64}}}

Reads a FORTRAN binary formatted abinit wavefunction file.

By default, abinit returns wavefunction data in `Complex{Float64}` entries, and this formatting is
maintained in the return type.

The header is used to automatically determine the file format, so this should read in any abinit
density output (provided a function exists to parse that header).
"""
read_abinit_wavefunction_new(filename) = open(read_abinit_wavefunction_new, filename)