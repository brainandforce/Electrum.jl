@testset "Atoms" begin
    a = NamedAtom("Cl1")
    @test a === NamedAtom("Cl1", 17)
    @test name(a) == "Cl1"
    @test atomic_number(a) === 17
    @test isdummy(a) === false
    @test NamedAtom(17) === NamedAtom("Cl", 17)
    @test isdummy(NamedAtom("test")) === true
end
