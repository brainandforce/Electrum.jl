# How to contribute to Electrum.jl

Thanks for taking interest in our package!

## Branch strategy

The `main` branch contains the latest development version of the code, which is guaranteed to build
and run (though it is possible it may not pass all tests). We also have a `next` branch which
contains code that has API-breaking changes which could affect code you've already written. The
current release version is kept at `release`, and older release versions can be accessed by
prefixing with the minor version number, like so: `0.1/release`.

As mentioned in the [README](README.md), you can track any of these branches with the Julia package
manager. In the future, you should be able to get Electrum.jl directly from the General registry.

### How to branch for contributions

If you'd like to implement a feature, create a new local branch, implement the feature, then push
branch and open a pull request on Github. To do this, you may need to fork the repo to your own
GitHub, and then push your changes there. The pull request tools on our instance of the repo should
let your create a pull request from your fork to our `main` or `next` branches.

We encourage you to prefix a branch with `feature/` if you are adding a feature to the package,
`fix/` for a fix, and `docs/` for changes to documentation. Note that features will only be merged
into `main`: old versions will not get new features, but they will get documentation updates and
fixes.

All merging should be done through a rebase operation. This is simpler than using a standard merge,
as it maintains a linear commit history. However, you will have to ensure that you have fetched or
pulled from the origin repo and rebased on the parent branch to avoid serious merge conflicts.

If you use VS Code, we recommend the 
[Git Graph](https://marketplace.visualstudio.com/items?itemName=mhutchie.git-graph) extension so you
can visualize what you're doing when you create a new branch, make commits, or push to a remote
repo.

## Continuous integration

Currently, we have automatically building documentation courtesy of Documenter.jl and its
integration with GitHub Actions. GitHub Actions are also used to automatically test the package for
every pull request and every merge into `main`. These tests should run automatically when you create
a pull request, and you should see the status of the tests after a few minutes. Pull requests that
do not pass tests generally will not be merged.

We use Aqua.jl to perform automatic package testing that covers some basic needs: method
ambiguities, stale dependencies, etc. Note that Aqua.jl is a bit picky about how `Project.toml` is
formatted.

## Dependencies

### Julia version

Electrum.jl is being written for Julia 1.6, which is an LTS release. This may change in the future,
but for now, avoid using any features that are present in later releases of Julia.

### Dependencies and interoperability

Avoid adding dependencies to Julia libraries that are not actively developed or maintained, or
contain functionality which would be simple to integrate into the package.

Try to minimize dependencies to code from different languages. Many other libraries pull functions
from scipy or other external libraries, but we intend to implement all core functionality here in
pure Julia.

## Coding style and standards

The following conventions are maintained throughout Electrum.jl.

When in doubt, follow the Julia style guide, located here:
https://docs.julialang.org/en/v1/manual/style-guide/

### Line length

Keep lines under 100 characters in all files. If a line is really long, consider whether you can
simplify the line by reducing the amount of nesting, shortening variable names, or splitting
operations into multiple lines.

For strings, we've found the best way (compatible with Julia 1.6+) is to call `string()` on a list
of strings to be concatentated.

### Naming

Stick to PascalCase for the names of types and modules, and snake_case for the names of variables
and functions, as with the majority of Julia code.

Functions for reading VASP inputs and outputs are usually given names like `readWAVECAR()`,
`writePOSCAR4()`,etc. with no spaces, and a version number afterwards for function writing. For
files which have the same format but a different name, for instance, `CONTCAR` files, you are
encouraged to add corresponding methods that contain the name of that particular file (for
instance, `readCONTCAR()`).

### File I/O functions

The [README.md](src/software/README.md) for the `src/software` directory provides notes on how to
write functions that perform file I/O.

### Comments and documentation

All functions and structs must have an associated docstring explaining the purpose of the code, 
*even if the struct or method is not exported.* For internal methods and structs, please prefix
them with the module name in the docstring.

It's better to be verbose about what's going on with your code. Even if a remark seems obvious, feel
free to leave it in. Perhaps the best assumption to make is that whoever is looking at the code may
not have any experience writing Julia code (or any code, for that matter).

Comments are generally placed above the lines that they refer to. Inline comments are fine; this is
just the pattern that's been used consistently in the code.

### Type parameters of newly defined types

All of the parametric types that contain dimensionality as a type parameter should have the
dimensionality parameters come first. So if you want to create `MyType` that has dimension `D` and
type `T`, the type should be created as `MyType{D,T}`. Dimension parameters are given as `D` in all
type definitions.

This is the opposite of the format used for Julia's built-in `AbstractArray` but matches the
convention used for `NTuple`.

### Constructors of parametric types

When defining a parametric type with an inner constructor for validation, ensure that the inner
constructor is a functor of a concrete type. For instance, for a struct `MyType{D,T}`, the 
declaration should look as follows:
```julia
struct MyType{D,T<:Number}
    f1::T
    f2::SVector{D,T}
    function MyType{D,T}(f1, f2::AbstractVector) where {D,T}
        # validation logic goes here
        return new(f1, f2)
    end
end
```
If some of the type parameters can be inferred from the data, create outer constructors that operate
on the partially parameterized types, like such.
```julia
function MyType{D}(f1, f2::AbstractVector) where D
    T = promote_type(typeof(f1), eltype(f2))
    return MyType{D,T}(f1, SVector{D}(f2))
end

function MyType(f1, f2::StaticVector)
    T = promote_type(typeof(f1), eltype(f2))
    return MyType{length(f2),T}(f1, f2)
end
```
### Type stability of constructors

Note that in the second outer constructor, `f2` must be a `StaticVector`. The reason why other
`AbstractVector` arguments are disallowed is because `length(f2)` can only be statically inferred
from a `StaticVector`. In general, only define type-stable constructors.

### Use of static vs. dynamic arrays

The `StaticArrays.jl` package provides support for static (fixed dimension) vectors, matrices, and
arrays. While dynamic arrays are convenient and don't require prior knowledge of array dimensions,
static arrays allow for higher performance.

The use of static arrays is encouraged whenever the array dimensionality is unlikely to change
(and in principle could be used as a type parameter). This universally applies to real space and
reciprocal space vectors, since their dimensionality is fixed (and usually 3).

The type `AtomList{D}` provides a great example of static vs. dynamic vector usage. The type
contains a `Vector{CartesianAtomPosition{D}}`, which is dynamic because the number of atoms in a
crystal (either the generating set or the visual template) may vary greatly. However, the type 
`AtomPosition{D}` contains a field for the atomic position vector, which is an `SVector{D,Float64}`
since atomic positions are not expected to vary in dimensionality.

The `StaticArrays` package provides the `SArray` type, which is immutable (cannot be altered after
creation), as well as the mutable `MArray` type. In general, it's better to stick to immutable 
types (though dynamic vectors are always mutable). For functions that should take only static array
arguments, use `::StaticArray{...}` in the type signature.

### `Tuple` and `NamedTuple`

When constructing a `Tuple` or a `NamedTuple`, parenthesize the construction for clarity. While
Julia can interpret statements like this correctly:
```julia
return a, b, c
```
it's a bit clearer to write it as
```julia
return (a, b, c)
```
Consider using a `NamedTuple` to work with related but disparate data. These can be thought of as
anonymous structs or immutable dictionaries that allow for indexing with a `Symbol` as well as a
numeric index:
```julia-repl
julia> nt = (a=1, b=2)
(a = 1, b = 2)

julia> nt[:b]
2
```
### Implementation of square matrix types

Julia currently does not allow field types to be computed. The problem with this is that `SMatrix`
must have a fourth type parameter that corresponds to the total length of the `NTuple` that is used
to construct it, even though that length can be inferred from the array dimensionality.

This poses a significant problem with the inclusion of matrices in struct defintions. If the type
assertion is `SMatrix{D1,D2,T}`, *the compiler will not be able to treat it as a concrete type.*
This means that small structs that could have otherwise been stack allocated end up pointing to data
on the heap, and this can reduce performance, especially in tight loops.

If you need to use an `SMatrix{D1,D2,T,L}` in a struct, be sure that you can define `L` as well as
`D1` and `D2`. Otherwise, it might be a better idea to store the data internally as an
`SVector{D2,SVector{D1,T}}`, and define `convert()` for it to turn it into a matrix. If you're
working with a mutable struct, this is not an issue: see `Electrum.Crystal` for an example of this
with the `transform` field. (`const` declarations in `mutable struct` cannot be used: they were
introduced in Julia 1.8, and we target backwards compatibility with Julia 1.6 LTS.)

### Logging and printing

Feel free to leave `@debug` statements in any code you include. Unless the logging level is set to
show them, they will not be seen by users. You can enable debug messages in the REPL for this
package by entering the following:

```julia-repl
julia> ENV["JULIA_DEBUG"] = Electrum
Electrum
```
Try to minimize the use of repeated `@info`, `@warn`, or `@error` statements in loops. Printing to
the terminal can bottleneck functions. `@debug` statements will normally be skipped unless a logger
explicitly compiles them or `$JULIA_DEBUG` is set to the module name. Alternatively, a function for
which explicit logging would be a useful feature can define a `quiet` keyword argument which can be
set to `true` for flexibility. Examples of this include `readWAVECAR` and `read_abinit_wfk`.

Avoid using `println()` if one of the logging macros better suits the purpose. If `println()` is
used, consider whether it might be a better idea to print to `stderr` instead of `stdout`. Note
that `@warn` and `@error` will print to `stderr` by default. (This means that if you try to pipe the
output to `less` or some other utility, it will not show up unless you perform output redirection.)
