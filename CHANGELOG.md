# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
  - New `EnergyOccupancy{T<:Real}` data structure for energy/occupancy pairs.
  - New `EnergiesOccupancies{T,N}` type alias for `Array{EnergyOccupancy{T},N}`. A constructor for
this type should be defined for all data structures containing energies and occupancies.
  - New functions: `energy`, `occupancy`, `energies`, `occupancies`
  - `EnergiesOccupancies` constructor for `PlanewaveWavefunction`

### Changed
  - `min_energy`, `max_energy`, `min_occupancy`, `max_occupancy`, and `fermi` all have generic
definitions based on `AbstractArray{<:EnergyOccupancy}`.
  - `PlanewaveWavefunction` uses these generic definitions via the constructor.

## [0.1.14]: 2023-11-02

### Added
  - New `RealBasis`/`ReciprocalBasis` constructors for `Electrum.ABINITHeader`.
  - New constructor for an empty `PlanewaveWavefunction` with known `KPointMesh`.

### Fixed
  - `read_abinit_WFK` now correctly stores the `KPointMesh` from the wavefunction file header.

## [0.1.13]: 2023-10-31

The developers would like to wish you a happy Halloween!

### Fixed
  - `readDOSCAR`, `readPROCAR`, and `get_fermi` now return energies in Hartree rather than eV.

## [0.1.12]: 2023-10-30

### Added
  - New `min_energy` and `max_energy` functions for `PlanewaveWavefunction`.
  - New internal `Electrum.dimension` function, which returns the type parameter of
`Electrum.DataSpace{D}` (usually a number). For objects and types with a `Electrum.DataSpace` trait,
this returns the associated dimension object.
  - Provided a fallback defintion for `Electrum.DataSpace(::Type{T}) where T`.

### Fixed
  - `BySpace{D}` did not subtype `DataSpace{D}`.

## [0.1.11]: 2023-10-22

### Added
  - New dependency on [InlineStrings.jl](https://github.com/JuliaStrings/InlineStrings.jl)

### Changed
  - `NamedAtom` now uses the `InlineString15` type instead of an `NTuple{8,Char}` for the `name`
field. Atom labels of up to 15 Unicode characters are supported.
  - `name(a::NamedAtom)` returns a `String`; `a.name` returns an `InlineStrings.InlineString15`, as
does `getfield(a, :name)`.

### Fixed
  - Admonition boxes in documentation now render correctly.

## [0.1.10]: 2023-10-19

### Added
  - Support for POSCAR and CONTCAR files with Cartesian coordinates.
  - All exported functions have docstrings, and those docstrings are included in the online
documentation.
  - `NamedAtom(::AbstractString)` now performs a more robust check for a valid atomic symbol. 
Strings such as `"Cl2"` are checked from the start to the first non-letter character.

### Changed
  - `readPOSCAR()` and `readCONTCAR()` now throw errors if the coordinate specification is missing.
  - Consolidated docstrings for several closely related types and functions, including
`readPOSCAR()` and `readCONTCAR()`, `Electrum.DataGrid` and its aliases, and `Electrum.LatticeBasis`
and its aliases.
  - `basis(x)` performs a type assertion to ensure an `Electrum.LatticeBasis` is returned.

## [0.1.9]: 2023-07-13

### Changed
  - Constructor for `KPoint{D}` does not automatically truncate the coefficients: this should
prevent bugs arising from imprecision in k-points with coefficients that are not exactly 0.5.

### Fixed
  - Reading of WAVECAR files should now index correctly: this was broken due to automatic rounding.

## [0.1.8]: 2023-07-13

### Fixed
  - Import of abinit binary outputs have the correct basis vectors associated with the data

## [0.1.7]: 2023-07-12

### Added
  - `Electrum.BySpace{D}` as a supertype for `ByRealSpace{D}` and `ByReciprocalSpace{D}`
  - `AbstractBasis{D,T}` as a type alias for `Electrum.LatticeBasis{<:Electrum.BySpace,D,T}`

### Changed
  - Cartesian indexing for `DataGrid` skips out on unnecessary bounds checking (due to periodicity).

## [0.1.6]: 2023-07-03

### Added
  - Support for reading abinit 9.10.1 binary outputs
  - Support for writing TOML files (optional thanks to Requires.jl dependency)

### Changed
  - `Electrum.get_abinit_version(::IO)` now strips spaces before constructing the version string

### Fixed
  - Some abinit docstrings had incorrect function names/arguments

## [0.1.5]: 2023-07-02

### Changed
  - When reading XSF files, the fractional coordinates are stored relative to the unit cell

### Fixed
  - Unit conversions weren't implemented for reading/writing XSF files

## [0.1.4]: 2023-06-28

The developers would like to wish you a happy [Tau Day](https://tauday.com/)!

### Added
  - Tests for XSF writing

### Fixed
  - More outstanding bugs for XSF writing

## [0.1.3]: 2023-06-27

### Fixed
  - Reference to old field name for previous data grid type was corrected

### Changed
 -  Added dependency on Requires.jl for optional features

## [0.1.2]: 2023-06-26

### Added
  - `partial_derivative()`, `gradient()`, and `directional_derivative()` functions for real space
data
  - `Electrum.SUnitVector{D,T}` data type for generating unit vectors. A unit vector of dimension
`D` and type `T` can be constructed with the unit at index `i` with `Electrum.SUnitVector{D,T}(i)`,
provided `oneunit(T)` is defined.

## [0.1.1] - 2023-06-21

### Changed
  - Unexported methods for reading anaddb inputs/outputs now have more standard signatures

### Fixed
  - `voxelsize(::DataGrid)` now correctly uses the grid length, not the basis vector length
  - `lengths(::LatticeBasis)` now correctly uses a column iterator to return the length values

## [0.1.0] - 2023-05-25

Initial release of Electrum.jl

[Unreleased]: https://github.com/brainandforce/Electrum.jl
[0.1.14]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.14
[0.1.13]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.13
[0.1.12]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.12
[0.1.11]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.11
[0.1.10]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.10
[0.1.9]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.9
[0.1.8]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.8
[0.1.7]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.7
[0.1.6]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.6
[0.1.5]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.5
[0.1.4]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.4
[0.1.3]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.3
[0.1.2]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.2
[0.1.1]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.0
