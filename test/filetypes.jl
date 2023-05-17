@testset "XSF files" begin
    # Check that the key name is correct
    @test haskey(xsf.data, "this_is_3Dgrid#1")
    # Check that the size of the RealSpaceDataGrid is correct
    @test size(xsf["this_is_3Dgrid#1"]) == (4, 4, 4)
end

@testset "abinit outputs" begin
    den = v80_den["density_total"]
    wfk = v80_wfk["wavefunction"]
    # Check that the correct FFT grid size is read
    @test size(den) == (24, 24, 36)
    # Check that there's 1 spin, 4 k-points, and 8 bands
    @test size(wfk)[1:3] == (1, 4, 8)
end

@testset "VASP outputs" begin
    @test length(poscar) == 4
    @test name.(atomtypes(poscar)) == ["Si", "Ir"]
    @test atomic_number.(atomtypes(poscar)) == [14, 77]
    @test natomtypes(poscar) == 2
    # Test that file writing works correctly
    writeCONTCAR(tmpdir, poscar)
    @test readCONTCAR(tmpdir) == sort(poscar)
end

@testset "LAMMPS position data" begin
    # This is going to suffer from floating point errors/parsing differences
    @test all(isapprox.(displacement.(lammps.atoms), displacement.(poscar.atoms), atol=1e-6))
    @test isapprox(basis(lammps), basis(poscar), atol=1e-6 * Electrum.ANG2BOHR)
end
