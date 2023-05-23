# Laves phase example
ScAl2_poscar = readPOSCAR("files/ScAl2.vasp")
# Transform from cF primitive -> conventional
F = [-1 1 1; 1 -1 1; 1 1 -1]

@testset "Crystals" begin
    ScAl2 = Crystal(ScAl2_poscar)
    @test convert(PeriodicAtomList, ScAl2) == ScAl2_poscar
    @test convert(PeriodicAtomList{3}, ScAl2) == ScAl2_poscar
    @test_throws TypeError convert(PeriodicAtomList{2}, ScAl2)
    @test basis(ScAl2) == basis(ScAl2_poscar)
    @test volume(basis(ScAl2)) == volume(basis(ScAl2_poscar))
    # Add transform for conventional cell
    set_transform!(ScAl2, F)
    @test ScAl2.transform == F
    @test generators(ScAl2) == ScAl2_poscar
    @test length(convert(PeriodicAtomList, ScAl2)) == length(ScAl2_poscar) * det(F)
    @test basis(ScAl2) == RealBasis{3}(basis(ScAl2_poscar))
    # Test equality
    ScAl2_cp = deepcopy(ScAl2)
    @test ScAl2 == ScAl2_cp
    @test ScAl2 !== ScAl2_cp                # different objects, same field values
    @test ScAl2 != Crystal(ScAl2_poscar)    # because the transform changed for ScAl2
end

@testset "Crystals with datasets" begin
    ScAl2_ds = CrystalWithDatasets(Crystal(ScAl2_poscar), Dict{String,Int}("test" => 1))
    @test propertynames(ScAl2_ds) === (:xtal, :data, :basis, :atoms, :sgno, :sgorig, :transform)
    set_transform!(ScAl2_ds, F)
    @test ScAl2_ds.xtal.transform == F
    @test ScAl2_ds.transform === ScAl2_ds.xtal.transform
    @test_throws TypeError convert(Crystal{2}, ScAl2_ds)
end
