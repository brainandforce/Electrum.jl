# Atoms

## Constructors and types
```@docs
NamedAtom
AbstractAtomPosition
CartesianAtomPosition
FractionalAtomPosition
AbstractAtomList
AtomList
PeriodicAtomList
```

## Extracting data
```@docs
name
atomic_number
isdummy
atomtypes
natomtypes
atomcounts
Base.isapprox(::T, ::T) where T<:AbstractAtomPosition
distance(::CartesianAtomPosition, ::CartesianAtomPosition)
```

## Moving and processing atom lists
```@docs
deduplicate
Electrum.move_into_cell
supercell
```
