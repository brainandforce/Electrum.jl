# Crystal data

## General
```@docs
Electrum.lengths(::Any)
Electrum.volume(::Any)
Electrum.basis(::Any)
```

## Basic data grids
```@docs
Electrum.AbstractDataGrid
Electrum.RealSpaceDataGrid
Electrum.voxelsize
Electrum.integrate
Electrum.HKLData
Electrum.fft
Electrum.ifft
```

## More complex data grids
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
