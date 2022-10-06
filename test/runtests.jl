using Test, Xtal

xsf = readXSF3D("files/test.xsf")

@testset "File formats" begin
    @test haskey("this_is_3Dgrid#1", xsf.data) broken=true
    @test size(xsf["this_is_3Dgrid#1           "]) == (4, 4, 4)
end
