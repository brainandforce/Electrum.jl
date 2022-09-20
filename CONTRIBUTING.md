# How to contribute to Xtal.jl

Thanks for taking interest in our package!

## Branch strategy

The `main` branch contains a version of the code that builds and can be used immediately. For now,
we have a `dev` branch that acts as a buffer while we test code manually. This will likely change
in the near future. We also have a `next` branch which contains code that has API-breaking changes
which could affect code you've already written. 

As mentioned in the [README](README.md), you can track any of these branches with the Julia package
manager if you'd like to test new features that have been added, or just feel like living life on
the edge.

If you'd like to implement a feature, create a new local branch, implement the feature, then push
branch and open a pull request on Github. To do this, you may need to fork the repo to your own
GitHub, and then push your changes there. The pull request tools on our instance of the repo should
let your create a pull request from your fork to our `main`, `dev`, or `next` branches.

All merging should be done through a rebase operation. This is simpler than using a standard merge,
as it maintains a linear commit history. However, you will have to ensure that you have fetched or
pulled from the origin repo and rebased on the parent branch to avoid serious merge conflicts.

## Testing and continuous integrations

As of now, there is not any testing implemented yet, nor are there continuous integrations. We hope
to use `Documenter` to generate documentation from docstrings in the near future.

## Dependencies

### Julia version

Xtal.jl is being written for Julia 1.6, which is an LTS release. This may change in the future, but
for now, avoid using any features that are present in later releases of Julia.

### Dependencies and interoperability

Avoid adding dependencies to Julia libraries that are not actively developed or maintained, or 
contain functionality which would be simple to integrate into the package.

Try to minimize dependencies to code from different languages. Many other libraries pull functions
from scipy or other external libraries, but we intend to implement all core functionality here in 
pure Julia.

## Coding style and standards

The following conventions are maintained throughout CrystalStructures.jl.

When in doubt, follow the Julia style guide, located here: https://docs.julialang.org/en/v1/manual/style-guide/

### Exceptions to the Julia style guide

For `Xtal.ABINITHeader` structs, direct field access may be used to assign and retrieve values.

Avoid using short-circuit `if` statements - write them out explicitly. There are many instances of 
short-circuited `if` statements throughout the code at this point, and you are welcome to restore 
them to full `if` statements.

### Line length

Keep lines under 100 characters in all files. You can use string concatenation, among other tools,
to split long strings if they come up. If a line is really long, consider whether you can simplify
the line by reducing the amount of nesting, shortening variable names, or splitting operations into
multiple lines. The general strategy that has been used is calling `string()`, which can stringify
a variety of arguments and pass through literals.

We have found that some methods of splitting strings may not be compatible with Julia 1.6 - please
make sure to test this before pushing.

### Naming

Stick to PascalCase for the names of types and modules, and snake_case for the names of variables
and functions.

### Comments and documentation

All functions and structs must have an associated docstring explaining the purpose of the code, 
*even if the struct or method is not exported.* For internal methods and structs, please prefix
them with the module name in the docstring.

It's better to be verbose about what's going on with your code. Even if a remark seems obvious, 
feel free to leave it in.

Comments are generally placed above the lines that they refer to. Inline comments are fine, this is
just the pattern that's been used consistently in the code.

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
reciprocal space vectors, since their dimensionality is fixed (and usually 3).

The type `AtomList{D}` provides a great example of static vs. dynamic vector usage. The type 
contains a `Vector{AtomPosition{D}}`, which is dynamic because the number of atoms in a crystal
(either the generating set or the visual template) may vary greatly. However, the type 
`AtomPosition{D}` contains a field for the atomic position vector, which is an `SVector{D,Float64}`
since atomic positions are not expected to vary in dimensionality.

The `StaticArrays` package provides the `SArray` type, which is immutable (cannot be altered after
creation), as well as the mutable `MArray` type. In general, it's better to stick to immutable 
types (though dynamic vectors are always mutable).

### Use of `BasisVectors{D}` and other reimplemented multidimensional array types

Julia currently does not allow field types to be computed. The problem with this is that `SMatrix`
must have a fourth type parameter that corresponds to the total length of the `NTuple` that is used
to construct it, even though that length can be inferred from the array dimensionality.

If you need to use an `SMatrix{D1,D2,T,L}` in a struct, be sure that you can define `L` as well as 
`D1` and `D2`. Otherwise, it might be a better idea to store the data internally as an
`SVector{D,SVector{D,T}}`, and define `convert()` for it to turn it into a matrix.

### Logging and printing

Feel free to leave `@debug` statements in any code you include. Unless the logging level is set to 
show them, they will not be seen by users. You can enable debug messages in the REPL for this
package by entering the following:

```julia-repl
julia> ENV["JULIA_DEBUG"] = Xtal
```

Try to minimize the use of repeated `@info`, `@warn`, or `@error` statements in loops. Printing to
the terminal can bottleneck functions. `@debug` statements will normally be skipped unless a logger
explicitly compiles them or `$JULIA_DEBUG` is set for the module.

Avoid using `println()` if one of the logging macros better suits the purpose. If `println()` is
used, consider whether it might be a better idea to print to `stderr` instead of `stdout`. Note 
that `@warn` and `@error` will print to `stderr` by default.
