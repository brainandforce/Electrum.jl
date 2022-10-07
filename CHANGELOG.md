# Changes from the 0.1 series

## New features

## Breaking changes

### Basis vectors and lattices

* `AbstractLattice` and its subtypes (`RealLattice` and `ReciprocalLattice`) have been removed.
* `BasisVectors` has been removed in favor of `RealBasis` and `ReciprocalBasis`, which are subtypes
of `AbstractBasis`. These have the advantage of explicitly declaring whether the vectors describe
real or reciprocal space.

### Crystal data
* The `pos` field of `Crystal`, which was used to store explicit atomic positions, has been 
removed, as have all methods that depend on this field.

### Pseudopotentials
* All potential and pseudopotential functionality has been removed from this package and will be 
placed in a new package.

### Renaming of structs, functions, and constants
* `cell_angle_cos()`, `cell_angle_rad()`, and `cell_angle_deg()` have been simplified to 
`angles_cos()`, `angles_rad()`, and `angles_deg()`.
* 
