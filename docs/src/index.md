# Electrum.jl

A Julia package that provides structs and methods for working with crystal structures.

# Introduction

Electrum.jl is a package designed to make the development of chemical theory tools easier. Not
only does it offer a type system that handles data commonly used in theory (such as real and
reciprocal space data grids), it also offers broad support for file types common to chemical theory.

In the future, we aim to support native Julia plotting of data processed by this package with a 
separate package that provides bindings to commonly used plotting utilities.

# Getting Electrum.jl

As of this release, Electrum.jl is not in the Julia package registry. You'll need to manually add
this repo to your Julia environment (which should be at least v1.6, and ideally the current release
version):

```
(@v1.8 pkg)> add https://github.com/brainandforce/Electrum.jl
```

You can also do this by importing `Pkg` and entering the following command:

```julia-repl
julia> import Pkg

julia> Pkg.add(url="https://github.com/brainandforce/Electrum.jl")
```

## Tracking different branches

If you'd like to work on developing Electrum.jl, you probably want to work on the project's current 
state, and not the release version. You can do this by specifying the branch you want to track
(let's assume it's `next`):

```
(@v1.8) pkg> add https://github.com/brainandforce/Electrum.jl#next
```
```julia-repl
julia> Pkg.add(url="https://github.com/brainandforce/Electrum.jl", rev="next")
```

# Licensing and attribution

Electrum.jl is [MIT licensed](https://mit-license.org/). In short, this means that you can use this
package however you like, without restrictions on relicensing.

You are not required to cite this package if you use it in research, however, attribution is
*always* appreciated. This is facilitated by the `CITATION.cff` file included in the repository.
