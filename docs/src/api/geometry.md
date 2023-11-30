# Geometric data

## Lattices

### Constructors and types
```@docs
Electrum.RealBasis
Electrum.ReciprocalBasis
Electrum.AbstractBasis
```

### Basis information from other data types
```@docs
Electrum.basis
Electrum.dualbasis
```

### Mathematical operations
```@docs
Base.inv(::Electrum.LatticeBasis)
Electrum.dual
Electrum.lengths(::Electrum.LatticeBasis)
Electrum.volume(::Electrum.LatticeBasis)
Electrum.angles_cos
Electrum.angles_rad
Electrum.angles_deg
Electrum.gram
Electrum.triangularize
```

### Iterators
```@docs
Electrum.eachvertex
```

### Miscellaneous
```@docs
Electrum.maxHKLindex
```

## Shift vectors
```@docs
Electrum.ShiftVector
Electrum.KPoint
Electrum.weight
```

## k-point meshes
```@docs
Electrum.KPointMesh
Electrum.nkpt
```
