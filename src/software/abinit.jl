#---Data structures--------------------------------------------------------------------------------#
"""
    Electrum.ABINITPseudopotentialInfo

Information about a pseudopotential used in an ABINIT calculation.
"""
struct ABINITPseudopotentialInfo
    title::String           # some info about the pseudopotential
    znuclpsp::Float64       # atomic number
    zionpsp::Float64        # charge
    pspso::Int32            # spin-orbit coupling considerations?
    pspdat::Int32           # revision date
    pspcod::Int32           # might be weird if 4 or 5, otherwise ignore
    pspxc::Int32            # XC functional used with pseudopotential
    lmn_size::Int32         # Not sure what this is 
    md5_pseudos::String      # Really, not as a UInt?
end

# Making this thing mutable is probably the best idea here
# That way we can initialize the struct and then update all the data
"""
    Electrum.ABINITHeader

Header information from an abinit FORTRAN binary output file. This data structure is mutable.
"""
Base.@kwdef mutable struct ABINITHeader
    # abinit version
    codvsn::VersionNumber = v"0.0.0"
    # Header format version
    headform::Int32 = 0
    # Unknown
    fform::Int32 = 0
    # Product of nkpt and nband
    bantot::Int32 = 0
    # Calculation date
    date::Int32 = 0
    # Grid interpolation scheme used for exchange-correlation energy                           
    intxc::Int32 = 0
    # Exchange-correlation functional
    ixc::Int32 = 0
    # Number of atoms
    natom::Int32 = 0
    # Size of the FFT grid
    ngfft::SVector{3,Int32} = zeros(SVector{3,Int32})
    # Number of k-points
    nkpt::Int32 = 0
    # Number of spin density components
    nspden::Int32 = 1
    # Number of spinor components
    nspinor::Int32 = 1
    # Number of spin polarizations
    nsppol::Int32 = 1
    # Number of symmetry operations
    nsym::Int32 = 1
    # Number of pseudopotentials used
    npsp::Int32 = 1
    # number of types of atoms
    ntypat::Int32 = 0
    # Method for determing band occupation
    occopt::Int32 = 0
    # Probably set to 1 if perturbative DFT is used
    pertcase::Int32 = 0
    # 1 if PAW is used
    usepaw::Int32 = 0                                           
    # Energy cutoff (in Hartree)
    ecut::Float64 = 0
    # Energy cutoff for the second grid (used in PAW)
    ecutdg::Float64 = 0                              
    # Smearing of the energy cutoff (used in optimizations)       
    ecutsm::Float64 = 0
    # ecut*dilatmx^2 (probably only set or different from ecut in optimizations)
    ecut_eff::Float64 = 0
    # q-points (needed for phonon calculations only)
    qptn::SVector{3,Float64} = zeros(SVector{3,Int32})
    # Translation vectors for the primitive cell
    # TODO: what exactly are the units for this thing?
    rprimd::SMatrix{3,3,Float64,9} = zeros(SMatrix{3,3,Float64})  
    # Bias voltage used in electron microscopy
    stmbias::Float64 = 0
    # Physical temperature of electrons (not sure of the units)
    tphysel::Float64 = 0
    # Smearing temperature for band occupation
    tsmear::Float64 = 0
    # Use of a wavelet basis set (we probably never use this)
    usewvl::Int32 = 0
    # Unknown for now
    nshiftk_orig::Int32 = 0                             # New for abinit 8.10
    # Number of k-point grid shifts
    nshiftk::Int32 = 0                                  # New for abinit 8.10
    # Maximum number of bands
    mband::Int32 = 0                                    # New for abinit 8.10
    # Set to 1 if the wavefunction should be stored at that k-point
    istwfk::Vector{Int32} = Int32[]    # size nkpt
    # Number of bands at each k-point
    nband::Vector{Int32} = Int32[]    # size nkpt
    # Number of planewave arrays at each k-point
    npwarr::Vector{Int32} = Int32[]    # size nkpt
    # Whether spin-orbit coupling is considered for each pseudopotential
    so_psp::Vector{Int32} = Int32[]    # size npsp
    # Antiferromagnetic symmetry components (-1 for those that reverse magnetization)
    symafm::Vector{Int32} = Int32[]    # size nsym
    # Symmetry operation matrices
    symrel::Vector{SMatrix{3,3,Int32,9}} = SMatrix{3,3,Int32}[]        # size nsym
    # Atom types (each integer corresponds to znucltypat)
    typat::Vector{Int32} = Int32[]    # size natom
    # List of k-points
    kpt::Vector{SVector{3,Float64}} = SVector{3,Float64}[]
    # Occupancies for each band in each k-point
    # TODO: changes this to a matrix or a vector of vectors
    occ::Vector{Float64} = Float64[]    # size bantot (nband*nkpt)
    # "Translation non-symmorphic vectors"
    tnons::Vector{SVector{3,Float64}} = SVector{3,Float64}[]    # size nsym
    # Atomic numbers for each nucleus
    znucltypat::Vector{Float64} = Float64[]    # size ntypat
    # k-point weights
    wtk::Vector{Float64} = Float64[]    # size nkpt
    # Setup method for k-point list
    kptopt::Int32 = 1
    # 2 if complex rho(i,j) occupancies are used (only for PAW calculations)
    pawcpxocc::Int32 = 1
    # Number of electrons in the cell
    nelect::Int32 = 0
    # Charge of the full cell
    cellcharge::Int32 = 0
    # 1 if wavelet solver is used/non-standard boundary conditions
    icoulomb::Int32 = 0
    # Lattice describing the k-point grid (ngkpt is along diagonal)
    kptrlatt::SMatrix{3,3,Int32,9} = zeros(SMatrix{3,3,Int32})
    # Unknown
    kptrlatt_orig::SMatrix{3,3,Int32,9} = zeros(SMatrix{3,3,Int32})
    # Unknown
    shiftk_orig::SVector{3,Float64} = zeros(SVector{3,Float64})
    # Shift of k-point grid off origin
    shiftk::SVector{3,Float64} = zeros(SVector{3,Float64})
    # Pseudopotential info
    pspinfo::Vector{ABINITPseudopotentialInfo} = ABINITPseudopotentialInfo[]    # size npsp
    # Residual density (a metric for calculation quality)
    residm::Float64 = 0
    # Positions of atoms (types given by typat matrix)
    xred::Vector{SVector{3,Float64}} = SVector{3,Float64}[]        # size natom
    # Total energy as determined by the calculation
    etotal::Float64 = 0
    # Fermi energy as determined by the calculation
    fermie::Float64 = 0
    # Masses of each type of atom (in amu)
    amu::Vector{Float64} = Float64[]
    # 2 if complex density information is used in PAW calculations
    cplex::Int32 = 1
    # rho(i,j) matrices for each atom
    pawdata::Matrix{Matrix{Float64}} = Matrix{Matrix{Float64}}(undef, 0, 0)    # size nspden*natom
end

# Equality definition
function Base.:(==)(h1::T, h2::T) where T<:ABINITHeader
    return all(getfield(h1, s) == getfield(h2, s) for s in fieldnames(T))
end

# Index notation for this thing, just in case that's easier
Base.getindex(h::ABINITHeader, name::Symbol) = getfield(h, name)
Base.setindex!(h::ABINITHeader, x, name::Symbol) = setfield!(h, name, x)

#---Constructors for commonly used types-----------------------------------------------------------#
RealBasis(h::ABINITHeader) = RealBasis(h.rprimd)
ReciprocalBasis(h::ABINITHeader) = ReciprocalBasis(RealBasis(h))

function Crystal(h::ABINITHeader)
    atomlist = PeriodicAtomList(
        RealBasis(h),
        FractionalAtomPosition.(Int.(h.znucltypat[h.typat]), h.xred)
    )
    # Don't include space group data since all atomic positions are generated (use defaults)
    return Crystal(atomlist)
end

function KPointMesh(h::ABINITHeader)
    (grid, shift) = (h.kptrlatt, h.shiftk)
    # abinit stores k-point weights normalized to 1 - convert this back to integers
    rweights = rationalize.(h.wtk / minimum(h.wtk); tol = sqrt(eps(Float64)))
    points = KPoint.(h.kpt, Int.(rweights .* maximum(x.den for x in rweights)))
    return KPointMesh(points, grid, shift)
end

#=
"""
    Electrum.symrel_to_sg(symrel::AbstractVector{<:AbstractMatrix{<:Integer}}) -> Int

Converts a list of ABINIT symmetry operations to the corresponding space group number and setting.
"""
function symrel_to_sg(symrel::AbstractVector{<:AbstractMatrix{<:Integer}})
    # any parity inverting operations?
    has_mir = !any(x -> x < 0, det.(symrel))
    # inversion center?
    has_inv = diagm([-1,-1,-1]) in symrel
    # All twofold axes
    twofold = [diagm([-1, -1, -1])]
    return 1 # for now
end
symrel_to_sg(h::ABINITHeader) = symrel_to_sg(h.symrel)
=#

# This function gets indices for the rhoij entries in a PAW calculation
function triang_index(n)
    row = floor(Int,(1 + sqrt(8*n - 7))/2)
    col = row + n - Int(row*(row + 1)/2)
    return (row, col)
end

#---Read header information------------------------------------------------------------------------#
"""
    Electrum.get_abinit_version(io::IO)
        -> NamedTuple{(:codvsn, :headform, :fform), Tuple{VersionNumber, Int32, Int32}}

Gets the ABINIT version information from a calculation output header.
"""
function get_abinit_version(io::IO)
    # Get the size of the first data entry
    sz = read(io, Int32)
    codvsn = VersionNumber(strip(String(read(io, sz - 8))))
    headform = read(io, Int32)
    fform = read(io, Int32)
    # Skip to the next field - assume the next thing to happen is one of the header reading methods
    skip(io, 4)
    return (codvsn=codvsn, headform=headform, fform=fform)
end

"""
    Electrum.read_abinit_header_57(io::IO) -> ABINITHeader

Reads in an abinit header from the outputs of calculations made by versions up to 7.10. These files
will contain a `headform` value of 57.

This function has been tested with outputs from abinit 7.10.5 using calculation data generated with
both norm-conserving pseudopotentials and PAW atomic data.
"""
function read_abinit_header_57(io::IO)
    # Variable names in this function correspond to names given here:
    # https://docs.abinit.org/guide/abinit/#header
    # Generate the new header
    h = ABINITHeader()
    # Skip the header for the first entry
    skip(io, 4)
    # Get a bunch of data now
    (
        h.bantot,     # total number of bands used in calculation (nband*nkpt*nsppol)
        h.date,       # date of calculation
        h.intxc,      # grid interpolation for XC
        h.ixc,        # XC functional index (usually 1 for Teter93 LDA)
        h.natom,      # number of atoms in unit cell
        ngfft1, ngfft2, ngfft3,    # FFT grid dimensions
        h.nkpt,       # number of k-points
        h.nspden,     # 1 for no spinors, 2 for scalar spin, 4 for vector spin
        h.nspinor,    # 1 for scalar wavefunction, 2 for spinor wavefunction
        h.nsppol,     # 1 for unpolarized, 2 for polarized
        h.nsym,       # number of symmetry operations
        h.npsp,       # number of pseudopotentials read
        h.ntypat,     # number of types of atoms
        h.occopt,     # occupation setting - usually 3
        h.pertcase,   # something about perturbations - probably not important
        h.usepaw      # 1 if the projector augmented wave method is used
    ) = tuple((read(io, Int32) for n = 1:18)...)
    # Reinterpret the FFT grid
    h.ngfft = SVector{3,Int}(ngfft1, ngfft2, ngfft3)
    # Get more data, but this time it's Float64s
    (
        h.ecut,       # Energy cutoff
        h.ecutdg,     # Energy cutoff for the second grid (used in PAW)
        h.ecutsm,     # Smearing of ecut for unit cell changes in optimization
        h.ecut_eff,   # ecut*dilatmx^2
    ) = tuple((read(io, Float64) for n in 1:4)...)
    # q-points (needed for phonon calculations only)
    h.qptn = SVector{3,Float64}((read(io, Float64) for n in 1:3)...)
    # primitive cell vectors in terms of acell
    # TODO: can we assume that the units are in Å for rprimd vectors?
    h.rprimd = SMatrix{3,3,Float64}((read(io, Float64) for n in 1:9)...)
    (
        h.stmbias,    # STM bias voltage: probably not used
        h.tphysel,    # Physical electron temperature
        h.tsmear      # Band smearing factor
    ) = tuple((read(io, Float64) for n in 1:3)...)
    h.usewvl = read(io, Int32)
    # next entry
    skip(io, 8)
    # storage of wavefunction by k-points
    h.istwfk = [read(io, Int32) for n in 1:h.nkpt]
    # number of bands (these elements should probably all be the same)
    h.nband = [read(io, Int32) for n in 1:(h.nkpt*h.nsppol)]
    # number of planewave arrays (or size ?)
    h.npwarr = [read(io, Int32) for n in 1:h.nkpt]
    # spin-orbit treatment of pseudopotentials (1 uses the psp file info)
    h.so_psp = [read(io, Int32) for n in 1:h.npsp]
    # symmetries for ferromagnetic materials (-1 if magnetization is altered)
    h.symafm = [read(io, Int32) for n in 1:h.nsym]
    # matrices for symmetry operations
    h.symrel = [SMatrix{3,3,Int32}(read(io, Int32) for n in 1:9) for n in 1:h.nsym]
    # types of atoms
    h.typat = [read(io, Int32) for n in 1:h.natom]
    # k-point coordinates
    h.kpt = [SVector{3,Float64}(read(io, Float64) for n in 1:3) for n in 1:h.nkpt]
    # Band occupancy
    h.occ = [read(io, Float64) for n in 1:h.bantot]
    # "translation non-symmorphic vectors"
    h.tnons = [SVector{3,Float64}(read(io, Float64) for n in 1:3) for n in 1:h.nsym]
    # Atomic numbers for each typat entry
    h.znucltypat = [read(io, Float64) for n in 1:h.ntypat]
    # k-point weights
    h.wtk = [read(io, Float64) for n in 1:h.nkpt]
    # Get all the pseudopotential info
    title = Vector{String}(undef, h.npsp)         # probably can be ignored
    znuclpsp = Vector{Float64}(undef, h.npsp)     # atomic number
    zionpsp = Vector{Float64}(undef, h.npsp)      # charge
    pspso = Vector{Int32}(undef, h.npsp)          # spin-orbit coupling considerations?
    pspdat = Vector{Int32}(undef, h.npsp)         # revision date
    pspcod = Vector{Int32}(undef, h.npsp)         # might be weird if 4 or 5, otherwise ignore
    pspxc = Vector{Int32}(undef, h.npsp)          # XC functional used with pseudopotential
                                                # probably matches ixc
    lmn_size = Vector{Int32}(undef, h.npsp)       # size of spherical harmonic components? Idk
    for m in 1:h.npsp
        skip(io, 8)
        title[m] = join(read(io, Char) for n in 1:132)
        (
            znuclpsp[m],
            zionpsp[m]
        ) = tuple((read(io, Float64) for n in 1:2)...)
        (
            pspso[m],
            pspdat[m],
            pspcod[m],
            pspxc[m],
            lmn_size[m]
        ) = tuple((read(io, Int32) for n in 1:5)...)
    end
    h.pspinfo = [
        ABINITPseudopotentialInfo(
            title[n],
            znuclpsp[n],
            zionpsp[n],
            pspso[n],
            pspdat[n],
            pspcod[n],
            pspxc[n],
            lmn_size[n],
            ""
        ) for n in 1:h.npsp
    ]
    # One more field
    skip(io, 8)
    # Residual density (gives the calculation quality)
    h.residm = read(io, Float64)
    # Atomic coordinates
    h.xred = [SVector{3,Float64}(read(io, Float64) for n in 1:3) for n in 1:h.natom]
    # Total energy and Fermi energy
    (h.etotal, h.fermie) = tuple((read(io, Float64) for n in 1:2)...)
    # Vector containing PAW data
    h.pawdata = [zeros(Float64, 8, 8) for s in 1:h.nspden, a in 1:h.natom]
    # Fields for a PAW calculation
    if h.usepaw == 1
        # First field
        skip(io, 8)
        # Number of entries in the rho_ij matrix (which is sparse and symmetric)
        nrhoijsel = [read(io, Int32) for s in 1:h.nspden, a in 1:h.natom]
        # Complex components of the density
        h.cplex = read(io, Int32)
        nsp_new = read(io, Int32)
        # Second field
        skip(io, 8)
        # Non-zero values of the rho_ij matrix
        rhoijselect = [
            [read(io, Int32) for n in 1:nrhoijsel[s,a]] for s in 1:h.nspden, a in 1:h.natom
        ]
        # Convert to matrix indices
        rhoijinds = [triang_index.(v) for v in rhoijselect]
        # Get the matrix elementwise
        # TODO: assuming that rhoijp consists of 8×8 matrices, but does it always?
        # rhoij = [zeros(Float64, 8, 8) for s in 1:nspden, a in 1:natom]
        for s in 1:h.nspden, a in 1:h.natom
            for ind in rhoijinds[s,a]
                x = read(io, Float64)
                h.pawdata[s,a][ind...] = x
                h.pawdata[s,a][reverse(ind)...] = x
            end
        end
    end
    # Skip the last record entry
    skip(io, 4)
    # Return output variables
    # TODO: perhaps try to calculate nelect assuming a neutral crystal?
    return h
end

"""
    Electrum.read_abinit_header_80(io::IO) -> ABINITHeader

Reads in an abinit header from the outputs of calculations made by versions up to 9.10. These files
will contain a `headform` value of 80.

This function has been tested with outputs from abinit 8.10.3 and abinit 9.10.1 using calculation
data generated with norm-conserving pseudopotentials.
"""
function read_abinit_header_80(io::IO)
    # Variable names in this function correspond to names given here:
    # https://docs.abinit.org/guide/abinit/#header
    # Generate the new header
    h = ABINITHeader()
    # Records start and end with the same Float32 containing the entry length
    skip(io, 4)                                           # seek to 26
    # Get a bunch of data now
    (
        h.bantot,     # total number of bands used in calculation (nband*nkpt*nsppol)
        h.date,       # date of calculation
        h.intxc,      # grid interpolation for XC
        h.ixc,        # XC functional index (usually 1 for Teter93 LDA)
        h.natom,      # number of atoms in unit cell
        ngfft1, ngfft2, ngfft3,    # FFT grid dimensions
        h.nkpt,       # number of k-points
        h.nspden,     # 1 for no spinors, 2 for scalar spin, 4 for vector spin
        h.nspinor,    # 1 for scalar wavefunction, 2 for spinor wavefunction
        h.nsppol,     # 1 for unpolarized, 2 for polarized
        h.nsym,       # number of symmetry operations
        h.npsp,       # number of pseudopotentials read
        h.ntypat,     # number of types of atoms
        h.occopt,     # occupation setting - usually 3
        h.pertcase,   # something about perturbations - probably not important
        h.usepaw      # 1 if the projector augmented wave method is used
    ) = tuple((read(io, Int32) for n = 1:18)...)
    # Reinterpret the FFT grid
    h.ngfft = SVector{3,Int}(ngfft1, ngfft2, ngfft3)
    # Get more data, but this time it's Float64s
    (
        h.ecut,       # Energy cutoff
        h.ecutdg,     # Energy cutoff for the second grid (used in PAW)
        h.ecutsm,     # Smearing of ecut for unit cell changes in optimization
        h.ecut_eff,   # ecut*dilatmx^2
    ) = tuple((read(io, Float64) for n in 1:4)...)
    # q-points (needed for phonon calculations only)
    h.qptn = SVector{3,Float64}((read(io, Float64) for n in 1:3)...)
    # primitive cell vectors in terms of acell
    # TODO: can we assume that the units are in Å for rprimd vectors?
    h.rprimd = SMatrix{3,3,Float64}((read(io, Float64) for n in 1:9)...)
    (
        h.stmbias,    # STM bias voltage: probably not used
        h.tphysel,    # Physical electron temperature
        h.tsmear      # Band smearing factor
    ) = tuple((read(io, Float64) for n in 1:3)...)
    (
        h.usewvl,
        h.nshiftk_orig,
        h.nshiftk,
        h.mband
    ) = tuple((read(io, Int32) for n in 1:4)...)
    # next entry
    skip(io, 8)
    # storage of wavefunction by k-points
    h.istwfk = [read(io, Int32) for n in 1:h.nkpt]
    # number of bands (these elements should probably all be the same)
    h.nband = [read(io, Int32) for n in 1:(h.nkpt*h.nsppol)]
    # number of planewave arrays (or size ?)
    h.npwarr = [read(io, Int32) for n in 1:h.nkpt]
    # spin-orbit treatment of pseudopotentials (1 uses the psp file info)
    h.so_psp = [read(io, Int32) for n in 1:h.npsp]
    # symmetries for ferromagnetic materials (-1 if magnetization is altered)
    h.symafm = [read(io, Int32) for n in 1:h.nsym]
    # matrices for symmetry operations
    h.symrel = [SMatrix{3,3,Int32}(read(io, Int32) for n in 1:9) for n in 1:h.nsym]
    # types of atoms
    h.typat = [read(io, Int32) for n in 1:h.natom]
    # k-point coordinates
    h.kpt = [SVector{3,Float64}(read(io, Float64) for n in 1:3) for n in 1:h.nkpt]
    # Band occupancy
    h.occ = [read(io, Float64) for n in 1:h.bantot]
    # "translation non-symmorphic vectors"
    h.tnons = [SVector{3,Float64}(read(io, Float64) for n in 1:3) for n in 1:h.nsym]
    # Atomic numbers for each typat entry
    h.znucltypat = [read(io, Float64) for n in 1:h.ntypat]
    # k-point weights
    h.wtk = [read(io, Float64) for n in 1:h.nkpt]
    # This is new to the version 80 header
    @debug string("position before new v80 data: ", position(io))
    # One more field
    skip(io, 8)
    # Residual density (gives the calculation quality)
    h.residm = read(io, Float64)
    # Atomic coordinates
    h.xred = [SVector{3,Float64}(read(io, Float64) for n in 1:3) for n in 1:h.natom]
    # Total energy and Fermi energy
    (h.etotal, h.fermie) = tuple((read(io, Float64) for n in 1:2)...)
    # Masses of atoms used
    h.amu = [read(io, Float64) for n in 1:h.ntypat]
    @debug string("position after amu: ", position(io))
    # This data is changed from version 57
    skip(io, 8)
    h.kptopt = read(io, Int32)
    h.pawcpxocc = read(io, Int32)
    h.nelect = read(io, Float64)
    h.cellcharge = read(io, Float64)
    h.icoulomb = read(io, Int32)
    # A lot of useful k-point info!
    @debug string("position before kptrlatt: ", position(io))
    h.kptrlatt = SMatrix{3,3,Int32}((read(io, Int32) for n in 1:9)...)
    h.kptrlatt_orig = SMatrix{3,3,Int32}((read(io, Int32) for n in 1:9)...)
    h.shiftk_orig = SVector{3,Float64}((read(io, Float64) for n in 1:3)...)
    h.shiftk = SVector{3,Float64}((read(io, Float64) for n in 1:3)...)
    # Get all the pseudopotential info
    title = Vector{String}(undef, h.npsp)           # probably can be ignored
    znuclpsp = Vector{Float64}(undef, h.npsp)       # atomic number
    zionpsp = Vector{Float64}(undef, h.npsp)        # charge
    pspso = Vector{Int32}(undef, h.npsp)            # spin-orbit coupling considerations?
    pspdat = Vector{Int32}(undef, h.npsp)           # revision date
    pspcod = Vector{Int32}(undef, h.npsp)           # might be weird if 4 or 5, otherwise ignore
    pspxc = Vector{Int32}(undef, h.npsp)            # XC functional used with pseudopotential
                                                    # probably matches ixc
    lmn_size = Vector{Int32}(undef, h.npsp)         # size of spherical harmonic components? Idk
    md5_pseudos = Vector{String}(undef, h.npsp)     # MD5 hash of pseudopotential used
    @debug string("before getting psp data: file pointer is at ", position(io))
    for m in 1:h.npsp
        skip(io, 8)
        @debug string(
            "Character in first field (should be H, probably):\n", 
            repr(MIME("text/plain"), peek(io, Char))
        )
        title[m] = join(read(io, Char) for n in 1:132)
        (
            znuclpsp[m],
            zionpsp[m]
        ) = tuple((read(io, Float64) for n in 1:2)...)
        (
            pspso[m],
            pspdat[m],
            pspcod[m],
            pspxc[m],
            lmn_size[m]
        ) = tuple((read(io, Int32) for n in 1:5)...)
        md5_pseudos[m] = String(read(io, 32))
    end
    @debug string("after getting psp data: file pointer is at ", position(io))
    h.pspinfo = [
        ABINITPseudopotentialInfo(
            title[n],
            znuclpsp[n],
            zionpsp[n],
            pspso[n],
            pspdat[n],
            pspcod[n],
            pspxc[n],
            lmn_size[n],
            md5_pseudos[n]
        ) for n in 1:h.npsp
    ]
    # Vector containing PAW data
    h.pawdata = [zeros(Float64, 8, 8) for s in 1:h.nspden, a in 1:h.natom]
    # Fields for a PAW calculation
    if h.usepaw == 1
        # First field
        skip(io, 8)
        # Number of entries in the rho_ij matrix (which is sparse and symmetric)
        nrhoijsel = [read(io, Int32) for s in 1:h.nspden, a in 1:h.natom]
        # Complex components of the density
        h.cplex = read(io, Int32)
        nsp_new = read(io, Int32)
        # Second field
        skip(io, 8)
        # Non-zero values of the rho_ij matrix
        rhoijselect = [
            [read(io, Int32) for n in 1:nrhoijsel[s,a]] for s in 1:h.nspden, a in 1:h.natom
        ]
        # Convert to matrix indices
        rhoijinds = [triang_index.(v) for v in rhoijselect]
        # Get the matrix elementwise
        # TODO: assuming that rhoijp consists of 8×8 matrices, but does it always?
        # rhoij = [zeros(Float64, 8, 8) for s in 1:nspden, a in 1:natom]
        for s in 1:h.nspden, a in 1:h.natom
            for ind in rhoijinds[s,a]
                x = read(io, Float64)
                h.pawdata[s,a][ind...] = x
                h.pawdata[s,a][reverse(ind)...] = x
            end
        end
    end
    # Skip the last record entry
    skip(io, 4)
    # Return output variables
    return h
end

"""
    Electrum.read_abinit_header(io::IO) -> Electrum.ABINITHeader

Reads the header of an ABINIT output file, automatically determining the format of the header from
the first few digits.
"""
function read_abinit_header(io::IO)
    # Get the info stored in the header
    (codvsn, headform, fform) = get_abinit_version(io)
    # @debug string("abinit version ", codvsn, " (header version ", headform, ")")
    # Read the rest of the header
    # Select the appropriate function by version number
    # header::ABINITHeader = fdict[headform](io)::ABINITHeader
    header = if headform == 57
        read_abinit_header_57(io)
    elseif headform == 80
        read_abinit_header_80(io)
    else
        error("Unsupported header format.")
    end
    # Add in the rest of the data and return
    header.codvsn = codvsn
    header.headform = headform
    header.fform = fform
    return header
end

read_abinit_header(filename) = open(read_abinit_header, filename)

#---Read datagrid information----------------------------------------------------------------------#
"""
    Electrum.read_abinit_datagrids(T, io, nspden, ngfft) -> Vector{Matrix{T}}

Reads the datagrid portion of an abinit density output, following the header.

Electron density values are given in electrons/Bohr³ in the abinit output files.

Depending on the value of `cplex` from the header, the resulting matrix may either contain
`Float64` or `Complex{Float64}` data.
"""
function read_abinit_datagrids(
    ::Type{T},
    io::IO,
    nspden::Integer,
    ngfft::AbstractVector{<:Integer};
    conversion::Real = 1
) where T
    # Density matrices
    rho = [Array{T,3}(undef, ngfft[1], ngfft[2], ngfft[3]) for s in 1:nspden]
    # Assuming that every grid is stored in its own 
    for s in 1:nspden
        @debug "Before assertion: file pointer is at " * string(position(io))
        sz1 = sizeof(T) * prod(ngfft)
        sz2 = read(io, Int32)
        @assert sz1 == sz2 string(
            "The reported FORTRAN formatted data entry doesn't match the grid size.\n",
            "Expected to read a datagrid of type ", T,
            " and size ", join(string.(ngfft), "×"), 
            " ( total size ", sz1 ,")\n", 
            "Reported size from file: ", sz2
        )
        @debug "After assertion: file pointer is at " * string(position(io))
        for a in 1:ngfft[1], b in 1:ngfft[2], c in 1:ngfft[3]
        # Convert to electrons/Ang^3
            rho[s][a,b,c] = read(io, T) * conversion
        end
        @debug "Read a datagrid: file pointer is at " * string(position(io))
        skip(io, 4)
    end
    return rho
end

#---Read whole files-------------------------------------------------------------------------------#
"""
    read_abinit_DEN(file)
        -> CrystalWithDatasets{3,String,RealDataGrid{3,Float64}}

Reads a FORTRAN binary formatted abinit density file. By default, abinit density files will have
the suffix `DEN`, but no assumptions are made about suffixes.

The header is used to automatically determine the file format, so this should read in any abinit
density output (provided a function exists to parse that header).

The number of datasets returned depends on the value of `nsppol` in the header, as calculations
with explicit treatment of spin will return spin densities.

Depending on the value of `cplex`, the datagrid(s) returned may be real or complex-valued.
"""
function read_abinit_DEN(io::IO)
    # Get the header from the file
    header = read_abinit_header(io)
    # Get the type of the datagrid (set by cplex in the header)
    T = (Float64, Complex{Float64})[header.cplex]
    rho = read_abinit_datagrids(T, io, header.nspden, header.ngfft)
    # Add each dataset to the dictionary
    data = Dict{String, RealDataGrid{3,T}}()
    # Convert the basis
    basis = RealBasis(header)
    # Fill the dictionary
    data["density_total"] = RealDataGrid(rho[1], basis)
    if header.nspden == 2
        data["density_spinup"] = RealDataGrid(rho[2], basis)
    elseif header.nspden == 4
        data["density_spinup_x"] = RealDataGrid(rho[2], basis)
        data["density_spinup_y"] = RealDataGrid(rho[3], basis)
        data["density_spinup_z"] = RealDataGrid(rho[4], basis)
    end
    return CrystalWithDatasets{3,String,RealDataGrid{3,T}}(Crystal(header), data)
end

read_abinit_DEN(filename) = open(read_abinit_DEN, filename)

"""
    read_abinit_POT(file)
        -> CrystalWithDatasets{3,String,RealDataGrid{3,T<:Number}}

Reads a FORTRAN binary formatted abinit potential file.

By default, abinit potential files will end in `POT` for the external potential, `VHA` for the 
Hartree potential, `VXC` for the exchange-correlation potential, and `VHXC` for the sum of both the
Hartree and exchange-correlation potentials.

The header is used to automatically determine the file format, so this should read in any abinit
density output (provided a function exists to parse that header).

The number of datasets returned depends on the value of `nsppol` in the header, as calculations with
explicit treatment of spin will return spin-dependent potentials.

Depending on the value of `cplex`, the datagrid(s) returned may be real or complex-valued.
"""
function read_abinit_POT(io::IO)
    # Get the header from the file
    header = read_abinit_header(io)
    # Get the type of the datagrid (set by cplex in the header)
    T = (Float64, Complex{Float64})[header.cplex]
    # No conversion will occur here: assume units of Hartree
    rho = read_abinit_datagrids(T, io, header.nspden, header.ngfft)
    data = Dict{String, RealDataGrid{3,T}}()
    basis = RealBasis(header)
    if header.nspden == 1
        data["potential_total"] = RealDataGrid(rho[1], basis)
    elseif header.nspden == 2
        data["potential_spinup"] = RealDataGrid(rho[1], basis)
        data["potential_spindown"] = RealDataGrid(rho[2], basis)
    elseif header.nspden == 4
        data["potential_up_up"] = RealDataGrid(rho[1], basis)
        data["potential_down_down"] = RealDataGrid(rho[2], basis)
        data["potential_up_down_real"] = RealDataGrid(rho[3], basis)
        data["potential_up_down_imag"] = RealDataGrid(rho[4], basis)
    end
    return CrystalWithDatasets{3,String,RealDataGrid{3,T}}(Crystal(header), data)
end

read_abinit_POT(filename) = open(read_abinit_POT, filename)

"""
    read_abinit_WFK(file; quiet = false)
        -> CrystalWithDatasets{3,String,PlanewaveWavefunction{3,Complex{Float64}}}

Reads a FORTRAN binary formatted abinit wavefunction file.

By default, abinit returns wavefunction data in `Complex{Float64}` entries, and this is maintained
in the return type.

The header is used to automatically determine the file format, so this should read in any abinit
density output (provided a function exists to parse that header).

By default, the function is verbose, with output printed for every k-point parsed, due to the large
size of the wavefunction files. If this behavior is undesirable, the `quiet` keyword argument may be
set to `true`.
"""
function read_abinit_WFK(io::IO; quiet = false)
    # Get the header from the file
    header = read_abinit_header(io)
    # Get the reciprocal lattice
    rlatt = ReciprocalBasis(header)
    # Get the minimum and maximum HKL values needed
    # Units for c (2*m_e/ħ^2) are hartree^-1 bohr^-2
    # this should affect only the size of the preallocated arrays
    bounds = SVector{3,UnitRange{Int}}(
        -g:g for g in maxHKLindex(rlatt, header.ecut)
    )
    @debug string(
        "hklbounds: ", bounds, "\n",
        "Calculated from ecut = ", header.ecut, " Hartree"
    )
    nb = maximum(header.nband)
    wf = PlanewaveWavefunction{3,Complex{Float64}}(
        rlatt,
        header.nsppol,
        KPointMesh(header),
        nb,
        bounds...
    )
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
            quiet || @info string(
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

read_abinit_WFK(filename; quiet = false) = open(f -> read_abinit_WFK(f; quiet), filename)

#---anaddb integrations (currently not part of the public API)-------------------------------------#

function read_abinit_anaddb_out(io::IO)
    # number of atoms
    readuntil(io, "natifc")
    num_atom = parse(Int64, readline(io))
    # fineness of grid
    # TODO: rename this variable: this doesn't seem to be a k-point mesh
    # But I'm not sure how to reconcile this with abinit ng2qpt
    readuntil(io, "ng2qpt")
    kptmesh = parse(Int64,split(readline(io))[1])
    # number of points along kpath
    readuntil(io,"nph1l")
    num_wavevect = parse(Int64,readline(io))
    # kpath - not used so commented out right now
    # kptlist = Vector{SVector{3,Float64}}(undef,num_wavevect)
    readuntil(io, "qph1l")
    readline(io)
    #==for wv in 1:num_wavevect
        kptlist[wv] = parse.(Float64,split(readline(io)))[1:3]
    ==#
    frequencies = zeros(Float64, num_wavevect, num_atom*3)
    # See later lines, but I think this should be a Matrix{SVector{3,Complex{Float64}}}
    modes = Matrix{SVector{6,Float64}}(undef, num_atom*3, num_atom)
    energies = Vector{Float64}(undef, num_atom*3)
    for wv in 1:num_wavevect
        readuntil(io, "in Thz")
        temp = split(readuntil(io, "Phonon energies"), ['\n',':',' '], keepempty=false)
        filter!(e->e≠"-",temp)
        frequencies[wv,:] = parse.(Float64, temp)
        # for now, only grab phonon modes at gamma, assumed to be the first point
        wv == 1 && for m in 1:num_atom*3
            readuntil(io, "Mode number")
            energies[m] = parse(Float64,split(readline(io))[3])
            for n in 1:num_atom
                temp = readline(io)
                # Skip lines if it warns about low frequency mode
                if startswith(strip(temp),"Attention")
                    temp = (readline(io); readline(io))
                end
                real = parse.(Float64,split(temp)[3:5])
                imag = parse.(Float64,split(readline(io))[2:4])
                # TODO: why not return an SVector{3,<:Complex}?
                # Perhaps the phonon modes could be their own type...
                modes[m,n] = SVector{6}(vcat(real, imag))
            end
        end
    end
    # TODO: can this use its own struct?
    return (
        #kptlist,
        kptmesh,
        FatBands{3}(
            frequencies,
            zeros(Float64, 9, num_atom, num_atom*3, num_wavevect),
            zeros(Complex{Float64},9, num_atom, num_atom*3, num_wavevect)
        ),
        modes,
        energies
    )
end

"""
    read_abinit_anaddb_out(file)
        -> Tuple{Int64, FatBands{3}, Array{SVector{6,Float64}}, Vector{Float64}}

Reads an output file from anaddb and returns the fineness of the kpoint path,
phonon dispersion bands (FatBands{3}), real and imaginary vectors for each atom
in each mode (Array{SVector{6,Float64}}) at gamma, and energies in Ha for each
mode (Vector{Float64}).

The fineness of the kpoint path is assumed to be the same for all three directions and
read from the first of the three numbers (Int64).

Phonon bands: Only the FatBands.bands field is filled with information. The other fields in FatBands
are defaulted to 0.

Modes: Each atom in the mode has real and imaginary vectors given in a static array, with the
first three numbers corresponding to the x, y, and z directions of the real vector, then
the next three numbers corresponding to the x, y, and z directions of the imaginary vector.
Note that the vectors are grabbed from the Gamma point in reciprocal space, which is assumed
to be the first k-point entered in the anaddb analysis.

Energies: Energies are given in a Vector{Float64} with a length matching the number of modes.
"""
read_abinit_anaddb_out(filename) = open(read_abinit_anaddb_out, filename)

function read_abinit_anaddb_in(io::IO)
    readuntil(io, "nph1l")  # number of phonons in list 1
    num_phonons = parse(Int, split(readline(io))[1])
    # num_phonons / kptmesh is used often, perhaps a variable binding is useful
    # Also where does kptmesh come from?
    kpath = Vector{String}(undef, ceil(num_phonons/kptmesh))
    readuntil(io, "qph1l")  # number of q-points in list 1
    for i in 1:floor(num_phonons/kptmesh)
        kpath[i] = filter(!isequal('!'), split(readline(io))[5])
        for j in 1:(kptmesh - 1)
            readline(io)
        end
    end
    # TODO: note why we need the 5th element
    kpath[ceil(num_phonons/kptmesh)] = filter(!isequal('!'), split(readline(io))[5])
    return kpath
end

"""
    read_abinit_anaddb_in(file) -> Vector{String}

Reads an input file for anaddb to determine the string of the path through reciprocal space
(`Vector{String}``) to plot the phonon dispersion curves. Also returns the fineness of the k-point
path.

Note that the strings we return gathered here are user-input. Example section of The
input file:
qph1l       0.00   0.00   0.00   1   !GAMMA
            0.00   0.05   0.00   1  
            0.00   0.10   0.00   1  
            0.00   0.15   0.00   1  
            0.00   0.20   0.00   1  
            0.00   0.25   0.00   1  
            0.00   0.30   0.00   1  
            0.00   0.35   0.00   1  
            0.00   0.40   0.00   1  
            0.00   0.45   0.00   1  
            0.00   0.50   0.00   1   !X

Here, the "GAMMA" and "X" are returned in the Vector{String}.
"""
read_abinit_anaddb_in(filename) = open(read_abinit_anaddb_in, filename)

# TODO: refactor this function into multiple methods
"""
    write_abinit_modes(
        modes::AbstractArray{<:StaticVector{6,<:Real}},
        energies::AbstractVector{<:Real}
    )

Writes the real and imaginary vectors into a .dat file for each mode. Each line in the .dat file
corresponds to the real x, y, z then imaginary x, y, z vectors of the corresponding atom. The
The vectors are Cartesian coordinates.
"""
function write_abinit_modes(
    modes::AbstractArray{<:StaticVector{6,<:Real}},
    energies::AbstractVector{<:Real}
)
    for i in axes(energies, 1)
        filename = string("mode_",i,"_energy_",energies[i],".dat")
        open(filename, "w") do io
            for n in axes(modes, 2)
                # prints real vector then imag vector
                # cartesian coordinates
                println(io, strip(string(modes[i,n]), ['[',']']))
            end
        end
    end
end

# TODO: rewrite this so that it doesn't allocate all of the lines
"""
    read_abinit_anaddb_PHDOS(file) -> Tuple{DensityOfStates, Vector{ProjectedDensityOfStates}}

Reads the PHDOS file from ABINIT's anaddb script and returns the total phonon density of states and
atomic contributions to the phonon density of states as a tuple. The PHDOS file is an output from
the anaddb program in ABINIT that contains information about the phonon density of states, with
frequencies in units of Hartree.
"""
function read_abinit_anaddb_PHDOS(file) 
    # Skip header (lines 1-8)
    data = readlines(file)[9:end]
    num_col = length(split(data[1]))
    # Parse data to floating point values
    # Energy | total PHDOS | int PHDOS | AtomType 1 PHDOS | AtomType 1 intPHDOS | ...
    data_new = [parse(Float64, split(data[i])[j]) for i in eachindex(data), j in 1:num_col]
    # Fermi for phonon density of states defaults to 0
    fermi = 0 
    tdos = DensityOfStates(fermi, data_new[:,1], data_new[:,2], data_new[:,3])
    pdos = [
        ProjectedDensityOfStates(
            fermi, 
            data_new[:,1],
            reshape(data_new[:,i*2+2], 1, size(data_new)[1])
        ) for i in 1:div((num_col-3), 2)
    ]
    return (tdos, pdos)
end
