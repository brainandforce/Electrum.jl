precompile(RealBasis, (AbstractMatrix{<:Real},))
precompile(ReciprocalBasis, (AbstractMatrix{<:Real},))

precompile(supercell, (PeriodicAtomList, Any))

precompile(readXSF, (AbstractString,))
precompile(writeXSF, (AbstractString, CrystalWithDatasets))

precompile(read_abinit_DEN, (AbstractString,))
precompile(read_abinit_POT, (AbstractString,))
precompile(read_abinit_WFK, (AbstractString,))

precompile(readPOSCAR, ())
precompile(readCONTCAR, ())
precompile(writePOSCAR, (Any,))
precompile(writeCONTCAR, (Any,))
precompile(readWAVECAR, ())
precompile(readDOSCAR, ())
