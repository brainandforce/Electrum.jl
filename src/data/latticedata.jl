"""
    AbstractGridPoints{D}

Supertype for data which specifies a set of grid points in `D`-dimensional space to be associated
with an array.
"""
abstract type AbstractGridPoints{D}
end

"""
    RegularGridPoints{D,T,B<:AbstractBasis{D,T},S<:StaticVector{D,T}}

Defines the placement of array data with respect to regularly space grid points within a lattice. 
This can be defined simply by wrapping an `AbstractBasis`, which assumes that data at the first
index of the array correspond to the origin point, and are regularly spaced along intervals equal to
the lengths of the lattice basis vectors divided by the number of grid points. Optional parameters
can be used to alter the straightforward assumption of placing data points:
  - The `shift` parameter moves the origin of the grid off of the origin coordinate of the lattice.
By default, this is the zero vector.
  - The `transform` parameter is a matrix of `Rational` numbers that describes how the data points
are transformed with respect to the basis vectors of the unit cell. By default, this is an identity
matrix, but other matrices may be useful when working with supercells of primitive cells.

When associating this type with an array, there are no assumptions in this data structure regarding
the size of the array.

# Type aliases

For convenience, the following aliases are defined:
    const RealRegularGridPoints{D,T} = RegularGridPoints{D,T,RealBasis{D,T},SVector{D,T}}
    const ReciprocalRegularGridPoints{D,T} = RegularGridPoints{D,T,ReciprocalBasis{D,T},KPoint{D,T}}

Note that in the case of `ReciprocalRegularGridPoints`, the shift vector is a `KPoint{D,T}`, whereas
for `RealRegularGridPoints` it is simply an `SVector{D,T}`. The reason for this is that data in
reciprocal space may be stored using weighted k-points
"""
struct RegularGridPoints{D,T,B<:AbstractBasis{D,T},S<:StaticVector{D,T}} <: AbstractGridPoints{D}
    basis::B
    shift::S
    _transform::SVector{D,{SVector{D,Rational{Int}}}}
    function RegularGridPoints{D,T,B,S}(
        basis,
        shift = zero(SVector{D,T}),
        transform = LinearAlgebra.I
    ) where {D,T,B}
        return new(basis, shift, SVector{D,Int}.(eachcol(transform)))
    end
end

const RealRegularGridPoints{D,T} = RegularGridPoints{D,T,RealBasis{D,T},SVector{D,T}}
const ReciprocalRegularGridPoints{D,T} = RegularGridPoints{D,T,ReciprocalBasis{D,T},KPoint{D,T}}

"""
    AbstractLatticeData{D,T,G<:AbstractGridPoints{D},A<:AbstractArray{T,D}} <: AbstractArray{T,D}

Supertype for arrays of type `A` containing data associated with spatial coordinates specified by
`G`.
"""
abstract type AbstractLatticeData{D,T,G<:AbstractGridPoints{D},A<:AbstractArray{T,D}} <: 
    AbstractArray{T,D}
end

"""
    Electrum.LatticeData{D,T,G<:RegularGridPoints,A} <: AbstractLatticeData{D,T,G,A}

Stores data in a `D`-dimensional array `A` with elements of type `T` defined on a regularly spaced
grid within a lattice `G`.

# Type aliases

Because in most instances, we know if we're working with real or reciprocal space data, we define
some aliases that refer specifically to those datasets:

    const RealLatticeData{D,T,X,A} = LatticeData{D,T,RealRegularGridPoints{D,X},A}
    const ReciprocalLatticeData{D,T,X,A} = LatticeData{D,T,ReciprocalRegularGridPoints{D,X},A}

The type parameter `X<:Real` is the shared element type of the lattice basis vectors and shift
vector.

We define the following aliases which are intended to be quick ways of referring to data wrapped by
a Julia `Array`:

    const RealLatticeArray{D,T,X} = RealLatticeData{D,T,X,Array{T,D}}
    const ReciprocalLatticeArray{D,T,X} = ReciprocalLatticeData{D,T,X,Array{T,D}}
"""
struct LatticeData{D,T,G<:RegularGridPoints,A} <: AbstractLatticeData{D,T,G,A}
    data::A
    grid::G
end

const RealLatticeData{D,T,X,A} = LatticeData{D,T,RealRegularGridPoints{D,X},A}
const ReciprocalLatticeData{D,T,X,A} = LatticeData{D,T,ReciprocalRegularGridPoints{D,X},A}

const RealLatticeArray{D,T,X} = RealLatticeData{D,T,X,Array{T,D}}
const ReciprocalLatticeArray{D,T,X} = ReciprocalLatticeData{D,T,X,Array{T,D}}

Base.size(l::LatticeData) = size(l.data)
Base.has_offset_axes(l::LatticeData) = true
Base.axes(l::LatticeData) = map(x -> x .- first(x), axes(l.data))

Base.getindex(l::LatticeData, i...) = @inbounds getindex(l.data, reinterpret_index(l, i)...)
Base.setindex!(l::LatticeData, x, i...) = @inbounds setindex!(l.data, x, reinterpret_index(l, i)...)

basis(l::LatticeData) = basis(l.grid)
