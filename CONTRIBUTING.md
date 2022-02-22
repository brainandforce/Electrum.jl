# How to contribute to Xtal.jl

Thanks for taking interest in our package! 

## Dependency versions

### Julia version

CrystalStructures.jl is being written for Julia 1.6, which is an LTS release. This may change in
the future, but for now, avoid using any features that are present in later releases of Julia.

## Coding standards

The following conventions are maintained throughout CrystalStructures.jl.

### Line length

Keep lines under 100 characters in all files. You can use string concatenation, among other tools,
to split long strings if they come up.

### Type parameters of newly defined types

All of the parametric types that contain dimensionality as a type parameter should have the 
dimensionality parameters come first. So if you want to create `MyType` that has dimension `D` and
type `T`, the type should be created as `MyType{D,T}`. Dimension parameters are given as `D` in all
type definitions.

This is the opposite of the format used for Julia's built-in `AbstractArray` but matches the 
convention used for `NTuple`.

### Use of static vs. dynamic arrays and tuples

The `StaticArrays.jl` package provides support for static (fixed dimension) vectors, matrices, and 
arrays. While dynamic arrays are convenient and don't require prior knowledge of array dimensions,
static arrays allow for higher performance.

The use of static arrays is encouraged whenever the array dimensionality is unlikely to change
(and in principle could be used as a type parameter). This universally applies to real space and
reciprocal space vectors, since their dimensionality is fixed.

The type `AtomList{D}` provides a great example of static vs. dynamic vector usage. The type 
contains a `Vector{AtomPosition{D}}`, which is dynamic because the number of atoms in a crystal
(either the generating set or the visual template) may vary greatly. However, the type 
`AtomPostion{D}` contains a field for the atomic position vector, which is an `SVector{D,Float64}`
since atomic positions are not expected to vary in dimensionality.

The `StaticArrays` package provides the `SArray` type, which is immutable (cannot be altered after
creation), as well as the mutable `MArray` type. In general, it's better to stick to immutable 
types (though dynamic vectors are always mutable).