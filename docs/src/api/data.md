# Crystal data

## General

## Basic data grids
```@docs
DataGrid
RealDataGrid
ReciprocalDataGrid
fft
ifft
AbstractFFTs.fftfreq(::DataGrid)
volume(::DataGrid)
voxelsize
integrate
partial_derivative
cell_gradient
gradient
directional_derivative
remove_shift 
```

## FFT iterators
```@docs
Electrum.FFTLength
Electrum.FFTBins
```

## Spin states
```@docs
Multiplicity
SpinBivector
```

## Energies and occupancies
```@docs
AbstractEnergyData
EnergyOccupancy
StateDensity
EnergiesOccupancies
energy
occupancy
density
energies
occupancies
densities
min_energy
max_energy
min_occupancy
max_occupancy
fermi
```

## Wavefunctions
```@docs
PlanewaveWavefunction
PlanewaveIndex
nspin
nband(::PlanewaveWavefunction)
nonzero_g_indices
nonzero_g_vectors
```

## Band structures
```@docs
BandAtKPoint
BandStructure
nband(::BandStructure)
nband(::BandAtKPoint)
FatBands
```

## Density of states
```@docs
AbstractDensityOfStates
DensityOfStates
ProjectedDensityOfStates
fermi(::AbstractDensityOfStates)
smear
energies(::AbstractDensityOfStates)
nelectrons
```

## Atomic data
```@docs
SphericalHarmonic
LinearAlgebra.dot(::SphericalHarmonic, ::SphericalHarmonic)
```
