# Geometric data

## Lattices

### Constructors and types
```@docs
Electrum.LatticeBasis
RealBasis
ReciprocalBasis
AbstractBasis
```

### Basis information from other data types
```@docs
basis
dualbasis
```

### Mathematical operations
```@docs
Base.inv(::Electrum.LatticeBasis)
dual
lengths(::Electrum.LatticeBasis)
volume(::Electrum.LatticeBasis)
angles_cos
angles_rad
angles_deg
gram
triangularize
```

### Iterators
```@docs
eachvertex
```

### Miscellaneous
```@docs
maxHKLindex
```

## Coordinate vectors
```@docs
CoordinateVector
RealCartesianCoordinate
RealFractionalCoordinate
ReciprocalCartesianCoordinate
ReciprocalFractionalCoordinate
ShiftVector
KPoint
weight
```

## k-point meshes
```@docs
KPointMesh
nkpt
```
