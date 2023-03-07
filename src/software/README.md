# Contribution notes for file handling

## Reading

When adding functions that read files, they need to be implemented in two methods. The first should
read from a file handle (Julia `IO` types) and return some object:
```julia
function read_file_format(io::IO; kwargs...)
    ...
end
```
The next one should apply generically to any other type of argument:
```julia
read_file_format(file; kwargs...) = open(io -> read_file_format(io; kwargs...), file)
```
Because Julia does not have a dedicated path type in `Base`, *string inputs should be interpreted as
file paths by default.* Path types provided by packages like
[FilePathsBase.jl](https://github.com/rofinn/FilePathsBase.jl) will overload `open()` with their own
path types, and the generic methods should work seamlessly with Electrum.

Try to keep the signature of the first method short (ideally, only accepting a file handle, with no 
keyword arguments). This keeps the implementation of the second method simpler.

## Writing

As above, file writing functions should consist of two methods.
```julia
function write_file_format(io::IO, args...; kwargs...)
    ...
end

function write_file_format(file, args...; kwargs...)
    open(file) do io
        write_file_format(io, args...; kwargs...)
    end
end
```
