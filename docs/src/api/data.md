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

## Wavefunctions
```@docs
Electrum.PlanewaveWavefunction
Electrum.PlanewaveIndex
Electrum.nspin
Electrum.nband(::PlanewaveWavefunction)
Electrum.fermi(::PlanewaveWavefunction)
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
Electrum.energies
Electrum.nelectrons
```

## Atomic data
```@docs
Electrum.SphericalHarmonic
LinearAlgebra.dot(::SphericalHarmonic, ::SphericalHarmonic)
```
