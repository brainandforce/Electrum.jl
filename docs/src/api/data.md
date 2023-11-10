# Crystal data

## General

## Basic data grids
```@docs
Electrum.DataGrid
Electrum.RealDataGrid
Electrum.ReciprocalDataGrid
Electrum.fft
Electrum.ifft
AbstractFFTs.fftfreq(::DataGrid)
Electrum.volume(::DataGrid)
Electrum.voxelsize
Electrum.integrate
Electrum.partial_derivative
Electrum.cell_gradient
Electrum.gradient
Electrum.directional_derivative
Electrum.remove_shift 
```

## FFT iterators
```@docs
Electrum.FFTLength
Electrum.FFTBins
```

## k-points
```@docs
Electrum.KPoint
Electrum.KPointMesh
Electrum.nkpt
Electrum.weight
```

## Energies and occupancies
```@docs
Electrum.EnergyOccupancy
Electrum.EnergiesOccupancies
Electrum.energy
Electrum.occupancy
Electrum.energies
Electrum.occupancies
Electrum.min_energy
Electrum.max_energy
Electrum.min_occupancy
Electrum.max_occupancy
Electrum.fermi
```

## Wavefunctions
```@docs
Electrum.PlanewaveWavefunction
Electrum.PlanewaveIndex
Electrum.nspin
Electrum.nband(::PlanewaveWavefunction)
Electrum.nonzero_g_indices
Electrum.nonzero_g_vectors
```

## Band structures
```@docs
Electrum.BandAtKPoint
Electrum.BandStructure
Electrum.nband(::BandStructure)
Electrum.nband(::BandAtKPoint)
Electrum.FatBands
```

## Density of states
```@docs
Electrum.AbstractDensityOfStates
Electrum.DensityOfStates
Electrum.ProjectedDensityOfStates
Electrum.fermi(::AbstractDensityOfStates)
Electrum.smear
Electrum.energies(::AbstractDensityOfStates)
Electrum.nelectrons
```

## Atomic data
```@docs
Electrum.SphericalHarmonic
LinearAlgebra.dot(::SphericalHarmonic, ::SphericalHarmonic)
```
