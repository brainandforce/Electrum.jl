# Lattices

The starting point for solid state structures is a periodic lattice. Electrum provides convenient
types for working with lattices of arbitrary dimension

## Real and reciprocal space traits

The `Electrum.BySpace{D}` supertype contains two types, `Electrum.ByRealSpace{D}` and
`Electrum.ByReciprocalSpace{D}`, where `D` is the number of dimensions associated with the dataset.
These types are used to denote whether data is associated with real space (e.g. electron density) or 
reciprocal space (e.g. the Fourier transform of the electron density). When working with lattices,
it is important to distinguish the two types of lattice: this is the primary reason why bare
`SMatrix{D,D,T}` instances are not used in this package.

## `Electrum.LatticeBasis` and methods

The `Electrum.LatticeBasis{S<:Electrum.BySpace,D,T}` data type is a wrapper for an
`SMatrix{D,D,T,D^2}` which represents the real or reciprocal space basis vectors of a lattice.

Electrum does not export `Electrum.LatticeBasis`, but instead provide the following aliases. This
allows developers to alter the implementation of `Electrum.LatticeBasis` without breaking the API:
```julia
const RealBasis = Electrum.LatticeBasis{ByRealSpace}
const ReciprocalBasis = Electrum.LatticeBasis{ByReciprocalSpace}
const AbstractBasis = Electrum.LatticeBasis{<:BySpace}
```
!!! warning
    Note that `RealBasis{D} === Electrum.LatticeBasis{ByRealSpace,D}`, and not
    `Electrum.LatticeBasis{ByRealSpace{D},D}`. While the latter is a valid type and can be used to
    store data, the result of `Electrum.DataSpace(::Electrum.LatticeBasis{ByRealSpace{D},D})` is
    an error!

The units of `RealBasis` are bohr, and those of `ReciprocalBasis` are radians over bohr,
corresponding with the convention that the dot product of a real basis vector with a corresponding
reciprocal basis vector is 2π.

### Construction

A `RealBasis` or `ReciprocalBasis` can be constructed from an `AbstractMatrix`, `Tuple`, or iterator
of the correct size.

In the case of a `StaticMatrix`, the size and element type are already known, so the constructors
can simply be called as `RealBasis` or `ReciprocalBasis`:
```
julia> RealBasis(SMatrix{3,3}(1, 0, 0, 0, 2, 0, 0, 0, 3))
Electrum.LatticeBasis{Electrum.ByRealSpace, 3, Int64}:
    a: [  1.000000   0.000000   0.000000 ]   (1.000000 bohr)
    b: [  0.000000   2.000000   0.000000 ]   (2.000000 bohr)
    c: [  0.000000   0.000000   3.000000 ]   (3.000000 bohr)
```
However, in the case of a `Matrix` or other dynamically sized data, the dimension is not known and
should be supplied to avoid an exception being thrown. This is to avoid type instability arising
from determining the size at runtime:
```
julia> RealBasis{3}([1 0 0; 0 2 0; 0 0 3])
Electrum.LatticeBasis{Electrum.ByRealSpace, 3, Int64}:
    a: [  1.000000   0.000000   0.000000 ]   (1.000000 bohr)
    b: [  0.000000   2.000000   0.000000 ]   (2.000000 bohr)
    c: [  0.000000   0.000000   3.000000 ]   (3.000000 bohr)
```
In either case, you can supply the element type (though it must be proceeded by the dimension in all
cases) to convert the input elements to the desired type (here, it's `Float32`):
```
julia> RealBasis{3,Float32}(SMatrix{3,3}(1, 0, 0, 0, 2, 0, 0, 0, 3))
Electrum.LatticeBasis{Electrum.ByRealSpace, 3, Float32}:
    a: [  1.000000   0.000000   0.000000 ]   (1.000000 bohr)
    b: [  0.000000   2.000000   0.000000 ]   (2.000000 bohr)
    c: [  0.000000   0.000000   3.000000 ]   (3.000000 bohr)
```
### Conversion

A `RealBasis` can be converted to a `ReciprocalBasis` via `Base.convert` or their constructors,
and vice versa:
```julia
b = ReciprocalBasis(a)
c = convert(RealBasis, b)
```
!!! tip
    Avoid needless back-and-forth conversions between `RealBasis` and `ReciprocalBasis` to avoid
    numerical instabilities.

### Mathematical operations

The `RealBasis` and `ReciprocalBasis` types support the majority of common operations used in 
solid-state chemistry, including addition, subtraction, multiplication, left division, and right
division. 

Importantly, it supports the QR decomposition provided by `LinearAlgebra.qr`. This decomposition is
useful in that it generates a `Q`` factor, which is an orthogonal matrix (representing a Euclidean
point isometry - compositions of rotations and reflections) and an ``R`` factor, which is an upper
triangular matrix. This operation is useful in converting lattices to a standard orientation: in the
case of a QR decomposition, the $R$ factor places the first basis vector of the lattice along the
first basis vector of space. In 3D, this means that ``\vec{a}`` is collinear with ``\vec{x}``.

Calling `LinearAlgebra.qr(::AbstractBasis{D,T})` returns a
`StaticArrays.QR{SMatrix{D,D,T,D^2}, SMatrix{D,D,T,D^2}, SVector{D,Int}}`, so the ``Q`` and ``R``
factors are bare ``SMatrix`` instances. While this is fine for the ``Q`` matrix, the ``R`` matrix
represents a basis, and we expect an `AbstractBasis` return value. The `triangularize` function
returns the ``R`` factor of a QR decomposition as an `AbstractBasis`, discarding the ``Q`` factor.

```@docs; canonical=false
Electrum.triangularize
```
!!! note
    Support for the corresponding LQ decomposition (where ``L`` is a lower triangular matrix) is not
    yet implemented.

### Implementation details

The `SMatrix{D1,D2,T,L}` type requires four type parameters - `D1` and `D2` are each of the matrix
dimensions, and `T` is the element type of the matrix. However, one last parameter is needed: `L`, 
the length of the `NTuple` that backs the `SMatrix`.

Julia currently does not allow for the calculation of type parameters from other type parameters,
which poses a problem in the declaration of structs. Technically, a type like `SMatrix{3,3,Float64}`
is an abstract type, as the `L` parameter is undeclared, even though the value of `L` is always
`D1 * D2`. By declaring a struct to have this type, there seems to be a performance drop.

The `LatticeBasis` types wrap an `SVector{D,SVector{D,Float64}}`. This only requires a single 
type parameter `D`, and allows for fully concrete struct declarations. The `:matrix` property
converts this data to an `SMatrix{D,D,T,D^2}` instance, and methods needing to access something that
looks like a matrix reference this property.

## Basis vectors in composite types

Some types (such as `Electrum.DataGrid`) store a set of basis vectors as part of the dataset. The
`basis` function allows a user to access that set of basis vectors. By default, `basis(x)` returns
`x.basis`, so defining a field `basis::RealBasis{...}` or `basis::ReciprocalBasis{...}` implements
these functions automatically.

!!! warning
    `basis(x)` returns an `AbstractBasis`, but there is no guarantee whether the return type is
    `RealBasis` or `ReciprocalBasis`. Use `convert(::Type{<:AbstractBasis}, basis(x))` to ensure
    that the return type is what you expect.

The type of basis vectors stored also allows for inference of the data space trait, as the default
definition of `DataSpace(::Type{T})` is `fieldtype(T, :basis)`. We encourage you to choose your
basis vector type to match the data you wish to represent.

!!! note
    For performance reasons, we encourage you to define struct fields with concrete types, and use 
    type parameters in your struct definitions to ensure that the type is concrete.
