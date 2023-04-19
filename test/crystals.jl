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
    @test volume(ScAl2_poscar) == volume(ScAl2)
    # Add transform for conventional cell
    set_transform!(ScAl2, F)
    @test ScAl2.transform == F
    @test generators(ScAl2) == ScAl2_poscar
    @test length(convert(PeriodicAtomList, ScAl2)) == length(ScAl2_poscar) * det(F)
    @test basis(ScAl2) == RealBasis{3}(basis(ScAl2_poscar) * F)
end

@testset "Crystals with datasets" begin
    ScAl2_dataset = CrystalWithDatasets(Crystal(ScAl2_poscar), Dict{String,Int}("test" => 1))
    @test propertynames(ScAl2_dataset) === (:xtal, :data, :atoms, :sgno, :sgorig, :transform,
        :slots, :keys, :vals, :ndel, :count, :age, :idxfloor, :maxprobe)
    set_transform!(ScAl2_dataset, F)
    @test ScAl2_dataset.xtal.transform == F
    @test ScAl2_dataset.transform === ScAl2_dataset.xtal.transform
    @test_throws TypeError convert(Crystal{2}, ScAl2_dataset)
end
