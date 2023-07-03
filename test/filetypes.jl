@testset "XSF files" begin
    # Check that the key name is correct
    @test haskey(xsf.data, "this_is_3Dgrid#1")
    g = xsf["this_is_3Dgrid#1"]
    # Check that the size of the RealSpaceDataGrid is correct
    @test size(g) == (4, 4, 4)
    # Check the indexing
    @test g[0] === 0.0
    @test g[1] === 1.000
    @test g[0,0,0] === 0.0
    @test g[1,1,1] === 1.732
    @test g[1,2,3] === 3.742
    # Test whether writing outputs and reading them back in gives us the same values
    output_filename = joinpath(tmpdir, "output.xsf")
    writeXSF(output_filename, xsf)
    xsf2 = readXSF(output_filename)
    h = xsf2["this_is_3Dgrid#1"]
    @test size(h) === size(g)
    @test iszero(h - g)
end

@testset "abinit outputs" begin
    den = v80_den["density_total"]
    wfk = v80_wfk["wavefunction"]
    # fform differs between the two - perhaps this could be useful later
    @test all(
        getfield(header_den, s) == getfield(header_wfk, s)
        for s in fieldnames(Electrum.ABINITHeader)[4:end]
    )
    @test header_den != header_wfk
    @test header_den[:fform] === header_den.fform
    @test Crystal(header_den) == Crystal(header_wfk)
    header_wfk[:fform] = header_den[:fform]
    @test header_den == header_wfk
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
    contcar = readCONTCAR(tmpdir)
    @test basis(contcar) ≈ basis(poscar)
    @test contcar.atoms == sort(poscar).atoms
end

@testset "LAMMPS position data" begin
    # This is going to suffer from floating point errors/parsing differences
    @test all(isapprox.(displacement.(lammps.atoms), displacement.(poscar.atoms), atol=1e-6))
    @test isapprox(basis(lammps), basis(poscar), atol=1e-6 * Electrum.ANG2BOHR)
end

@testset "TOML" begin
    using TOML
    xtal = v80_den.xtal
    toml_path = joinpath(tmpdir, "test.toml")
    writeTOML(toml_path, xtal)
    toml = TOML.parsefile(joinpath(tmpdir, "test.toml"))
    @test toml["transform"] == collect(eachcol(xtal.transform))
    @test toml["basis"]["vectors"] == collect(eachcol(basis(xtal)))
    @test toml["basis"]["dimension"] == 3
    @test toml["basis"]["realspace"]
    @test length(toml["atoms"]) == 2
    @test toml["atoms"][1]["name"] == "Sc"
end
