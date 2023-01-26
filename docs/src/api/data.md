# Crystal data

## General
```@docs
Electrum.lengths(::Any)
Electrum.volume(::Any)
```

## Real space
```@docs
Electrum.AbstractRealSpaceData
Electrum.RealSpaceDataGrid
Electrum.grid
Electrum.voxelsize
Electrum.integrate
Electrum.fft
```

## Reciprocal space
```@docs
Electrum.AbstractReciprocalSpaceData
Electrum.AbstractKPoints
Electrum.AbstractHKL
Electrum.KPointGrid
Electrum.KPointList
Electrum.BandAtKPoint
Electrum.BandStructure
Electrum.FatBands
Electrum.HKLData
Electrum.HKLDict
Electrum.ReciprocalWavefunction
Electrum.nkpt
Electrum.nband
Electrum.ifft
```

## Density of states
```@docs
Electrum.AbstractDensityOfStates
Electrum.DensityOfStates
Electrum.ProjectedDensityOfStates
Electrum.nkpt
Electrum.nband
Electrum.bounds
Electrum.fermi
Electrum.smear
Electrum.energies
Electrum.nelectrons
```

## Atomic data
```@docs
Electrum.SphericalHarmonic
```
