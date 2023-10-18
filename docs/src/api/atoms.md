# Atoms

## Constructors and types
```@docs
Electrum.NamedAtom
Electrum.AbstractAtomPosition
Electrum.CartesianAtomPosition
Electrum.FractionalAtomPosition
Electrum.AbstractAtomList
Electrum.AtomList
Electrum.PeriodicAtomList
```

## Extracting data
```@docs
Electrum.name
Electrum.atomic_number
Electrum.isdummy
Electrum.atomtypes
Electrum.natomtypes
Electrum.atomcounts
Base.isapprox(::T, ::T) where T<:AbstractAtomPosition
```

## Moving and processing atom lists
```@docs
Electrum.deduplicate
Electrum.move_into_cell
Electrum.supercell
```
