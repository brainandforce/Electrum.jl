# Changes from the 0.1 series

## New features

## Breaking changes

* All potential and pseudopotential functionality has been removed from this package and will be 
placed in a new package, `DFTPotentials.jl`.
* Many unnecessarily long function names have been removed.
** `cell_angle_cos()`, `cell_angle_rad()`, and `cell_angle_deg()` have been simplified to 
`angles_cos()`, `angles_rad()`, and `angles_deg()`.
** 

## Organizational changes

The `data.jl` file will be split into multiple files within the `src/data` directory. (This has
caused some of the file history to be lost to Git)
