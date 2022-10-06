# Xtal.jl

A Julia package that provides structs and methods for working with crystal structures.

# Introduction

Xtal.jl is a package designed to make the development of chemical theory tools easier. Not only 
does it offer a type system that handles data commonly used in theory (such as real and reciprocal
space data grids), it also offers broad support for filetypes common to chemical theory.

In the future, we aim to support native Julia plotting of data processed by this package with a 
separate package.

# Getting Xtal.jl

As of this release, Xtal.jl is not in the Julia package registry. You'll need to manually add this
repo to your Julia environment (which should be at least v1.6):

```
(@v1.8 pkg)> add https://github.com/brainandforce/Xtal.jl
```

You can also do this by importing `Pkg` and entering the following command:

```julia-repl
julia> import Pkg

julia> Pkg.add(url="https://github.com/brainandforce/Xtal.jl")
```

## Tracking different branches

If you'd like to work on developing Xtal.jl, you probably want to work on the project's current 
state, and not the release version. You can do this by specifying the branch you want to track
(let's assume it's `dev`):

```
(@v1.8) pkg> add https://github.com/brainandforce/Xtal.jl#dev
```
```julia-repl
julia> Pkg.add(url="https://github.com/brainandforce/Xtal.jl", rev="dev")
```

# License

Xtal.jl is MIT licensed.
