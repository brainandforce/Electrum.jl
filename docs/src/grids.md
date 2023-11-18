# Data grids

The output of many quantum chemistry packages are datasets defined on a spatial grid in real or
reciprocal space. Theoretical tools may want to perform mathematical operations on this data, but
operating on bare arrays containing the dta requires tracking properties associated with each array,
such as the basis vectors of the lattice

## `DataGrid`, `RealDataGrid`, and `ReciprocalDataGrid`

An `DataGrid{D,B<:Electrum.LatticeBasis,S<:AbstractVector{<:Real},T}` contains data defined in a
crystal lattice of `D` with basis vectors of type `B`, a shift parameter of type `S`, and elements
of type `T`, either in real space (`RealDataGrid`) or in reciprocal space (`ReciprocalDataGrid`).
These aliases are defined as follows:

```julia
const RealDataGrid{D,T} = DataGrid{D,RealBasis{D,Float64},SVector{D,Float64},T}
const ReciprocalDataGrid{D,T} = DataGrid{D,ReciprocalBasis{D,Float64},KPoint{D},T}
```
!!! note
    Look closely at the above definition: the shift data type for `RealDataGrid{D}` is
    `SVector{D,Float64}`, but the shift data type for `ReciprocalDataGrid{D}` is `KPoint{D}`. When
    the `DataGrid` constructor without the shift type parameter is invoked, the basis type is used
    to infer the appropriate shift type so that a `RealDataGrid` or `ReciprocalDataGrid` is
    constructed.

`DataGrid` uses zero-based, periodic indexing: the first index of an `AbstractDataGrid{D}` is
`zero(NTuple{D,Int})`, and indices whose moduli with respect to size along that dimension are
identical will reference the same element: for instance, for `g::AbstractDataGrid{3}` with size
`(10, 10, 10)`, `g[69, 420, 1337] === g[9, 0, 7]`. Encountering a `BoundsError` is not possible
when indexing an `DataGrid`.

The basis of an `DataGrid` can be recovered with `basis(::DataGrid)`, which will be of the type
specified by the type parameter. 

## Broadcasting and mathematical operations

Broadcasting is defined for `DataGrid` with a custom `Base.Broadcast.BroadcastStyle` subtype:
```julia
Electrum.DataGridStyle{D,B,S} <: Broadcast.AbstractArrayStyle{D}
```
This allows `DataGrid` instances to operated on with dot syntax. However, they must share lattice
basis vectors and shift values. If they do not match, an `Electrum.LatticeMismatch` exception will
be thrown.

!!! info
    Although `Base.Broadcast.ArrayStyle` is usually overridden by other subtypes of 
    `Base.Broadcast.AbstractArrayStyle`, it does not override `Electrum.DataGridStyle`. Adding a
    `DataGrid` to an `Array` returns an `Array`, and adding a `DataGrid` to other `AbstractArray`
    subtypes returns the `AbstractArray` subtype defined by the `Broadcast.BroadcastStyle`. In the 
    case of a dimension mismatch, the broadcast style wll be `Broadcast.ArrayConflict` - the
    operation will throw a `DimensionMismatch`.

The `+` and `-` operators are defined for `DataGrid` instances, and they are faster than the
broadcasted `.+` and `.-` equivalents. As with the broadcasted versions, checks are implemented to
ensure that the lattice basis vectors and shifts match.

Similarly, the `*`, `/`, and `\` operators are defined for pairs of `DataGrid` and `Number`
instances, and again, are faster than their broadcasted equivalents.

### Fourier transforms

The Fourier transform and its inverse are available through an overload of `FFTW.fft()` and
`FFTW.ifft()`. The transforms are normalized with respect to the basis vectors of the space, so
for `g::DataGrid`, `ifft(fft(g)) ≈ g` (to within floating point error).

`fft` and `ifft` defined for `DataGrid` generally, but it is critical to note that in the vast
majority of cases, you will want to call `fft` on `RealDataGrid` instances and `ifft` on
`ReciprocalDataGrid` instances.
