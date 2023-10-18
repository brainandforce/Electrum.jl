# Crystal data

## General

## Basic data grids
```@docs
Electrum.DataGrid
Electrum.RealDataGrid
Electrum.ReciprocalDataGrid
Electrum.fft
Electrum.ifft
Electrum.volume(::DataGrid)
Electrum.voxelsize
Electrum.integrate
Electrum.partial_derivative
Electrum.gradient
Electrum.directional_derivative
Electrum.remove_shift 
```

## Wavefunctions
```@docs
Electrum.PlanewaveWavefunction
Electrum.PlanewaveIndex
Electrum.nspin(::PlanewaveWavefunction)
Electrum.nkpt(::PlanewaveWavefunction)
Electrum.nband(::PlanewaveWavefunction)
Electrum.fermi(::PlanewaveWavefunction)
```

## k-points
```@docs
Electrum.KPoint
Electrum.KPointMesh
```

## Band structures
```@docs
Electrum.BandAtKPoint
Electrum.BandStructure
Electrum.nkpt(::BandStructure)
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
Electrum.energies
Electrum.nelectrons
```

## Atomic data
```@docs
Electrum.SphericalHarmonic
```
