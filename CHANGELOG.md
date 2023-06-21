# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
 -  Unexported methods for reading anaddb inputs/outputs now have more standard signatures

### Fixed
 -  `voxelsize(::DataGrid)` now correctly uses the grid length, not the basis vector length
 -  `lengths(::LatticeBasis)` now correctly uses a column iterator to return the length values

## [0.1.0] - 2023-05-25

Initial release of Electrum.jl

[Unreleased]: https://github.com/brainandforce/Electrum.jl
[0.1.1]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/brainandforce/Electrum.jl/releases/tag/v0.1.0
