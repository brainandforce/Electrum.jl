using Test, Xtal

xsf = readXSF3D("files/test.xsf")

@testset "File formats" begin
    @test haskey(xsf.data, "this_is_3Dgrid#1")
    @test size(xsf["this_is_3Dgrid#1"]) == (4, 4, 4)
end
