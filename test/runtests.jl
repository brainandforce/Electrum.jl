using LinearAlgebra, StaticArrays
using Test, Aqua, Xtal

Aqua.test_all(Xtal; project_toml_formatting=false)

xsf = readXSF3D("files/test.xsf")
v80_den = read_abinit_density("files/Sc_eq_o_DEN")
v80_wfk = read_abinit_wavefunction("files/Sc_eq_o_WFK")
poscar = readPOSCAR("files/")

@testset "All tests" begin
    @testset "Basis vectors" begin
        # Get basis vectors from the XSF file
        b = basis(xsf["this_is_3Dgrid#1"])
        # Inversion should give us 2π along the diagonal
        @test ReciprocalBasis(b) == ReciprocalBasis{3}(diagm([2π,2π,2π]))
        # Check that conversion between real and reciprocal bases is invertible
        @test b ≈ RealBasis(ReciprocalBasis(b))
    end

    @testset "Supercell construction" begin
        # Coordinates for Laves-MgZn2
        MgZn2 = AtomList(
            RealBasis(5.223 * SMatrix{3,3,Float64}(1, 0, 0, -1/2, sqrt(3/4), 0, 0, 0, sqrt(8/3))),
            [
                AtomPosition(12, SVector(1/3, 2/3, 1/16)),
                AtomPosition(12, SVector(1/3, 2/3, 7/16)),
                AtomPosition(12, SVector(2/3, 1/3, 9/16)),
                AtomPosition(12, SVector(2/3, 1/3, 15/16)),
                AtomPosition(30, SVector(0.0, 0.0, 0.0)),
                AtomPosition(30, SVector(0.0, 0.0, 1/2)),
                AtomPosition(30, SVector(5/6, 2/3, 1/4)),
                AtomPosition(30, SVector(5/6, 1/6, 1/4)),
                AtomPosition(30, SVector(1/3, 1/6, 1/4)),
                AtomPosition(30, SVector(1/6, 1/3, 3/4)),
                AtomPosition(30, SVector(1/6, 5/6, 3/4)),
                AtomPosition(30, SVector(2/3, 5/6, 3/4)),
                # Two duplicates for good measure!
                AtomPosition(12, SVector(4/3, 2/3, 1/16)),
                AtomPosition(30, SVector(0.0, 0.0, 1/2))
            ]
        )
        # The extra Zn atom should be stripped because it's identical to an earlier one
        @test length(MgZn2) == 13
        # Build an orthorhombic supercell (treat hexagonal as orthorhombic C)
        transform = [1 1 0; -1 1 0; 0 0 1]
        MgZn2_sc = supercell(MgZn2, transform)
        # Check that the basis is actually orthorhombic
        @test convert(SMatrix, basis(MgZn2_sc)) ≈ diagm(SVector(9.0465013679, 5.223, 8.5291232843))
        # Check that the atomic list size scales with the transform matrix
        @test length(MgZn2_sc)/length(MgZn2) <= abs(det(transform))
        # Check that the duplicate Mg atom was stripped
        @test length(MgZn2_sc) == 24
        # Check each individual atom:
        #=
        @testset "MgZn2 supercell construction check" begin
            # Just work with the coordinates
            Mg = [coord(a) for a in filter(12, MgZn2_sc)]
            Zn = [coord(a) for a in filter(30, MgZn2_sc)]
            @test any(isapprox(SVector(1/2,   1/6,  1/16)), Mg)
            @test any(isapprox(SVector(1/2,   1/6,  7/16)), Mg)
            @test any(isapprox(SVector(1/2,   5/6,  9/16)), Mg)
            @test any(isapprox(SVector(1/2,   5/6, 15/16)), Mg)
            @test any(isapprox(SVector(  0,   2/3,  1/16)), Mg)
            @test any(isapprox(SVector(  0,   2/3,  7/16)), Mg)
            @test any(isapprox(SVector(  0,   1/3,  9/16)), Mg)
            @test any(isapprox(SVector(  0,   1/3, 15/16)), Mg)
            @test any(isapprox(SVector(0.0,   0.0,   0.0)), Zn)
            @test any(isapprox(SVector(  0,     0,   1/2)), Zn)
            @test any(isapprox(SVector(1/2,   1/2,   0.0)), Zn)
            @test any(isapprox(SVector(1/2,   1/2,   1/2)), Zn)
            @test any(isapprox(SVector(  0,   1/6,   1/4)), Zn)
            @test any(isapprox(SVector(  0,   5/6,   3/4)), Zn)
            @test any(isapprox(SVector(1/4,  5/12,   1/4)), Zn)
            @test any(isapprox(SVector(3/4,  5/12,   1/4)), Zn)
            @test any(isapprox(SVector(1/2,   1/3,   1/4)), Zn)
            @test any(isapprox(SVector(1/4, 11/12,   1/4)), Zn)
            @test any(isapprox(SVector(3/4, 11/12,   1/4)), Zn)
            @test any(isapprox(SVector(1/4,  1/12,   3/4)), Zn)
            @test any(isapprox(SVector(3/4,  1/12,   3/4)), Zn)
            @test any(isapprox(SVector(1/2,   1/3,   3/4)), Zn)
            @test any(isapprox(SVector(1/4,  7/12,   3/4)), Zn)
            @test any(isapprox(SVector(3/4,  7/12,   3/4)), Zn)
        end
        =#
    end

    @testset "RealSpaceDataGrid" begin
        g = xsf["this_is_3Dgrid#1"]
        # Check the indexing
        @test g[0] === 0.0
        @test g[1] === 1.000
        @test g[0,0,0] === 0.0
        @test g[1,1,1] === 1.732
        @test g[1,2,3] === 3.742
    end

    @testset "Fourier transforms" begin
        g = xsf["this_is_3Dgrid#1"]
        hkl = fft(g)
        # Check that the Fourier transform and its inverse undo each other
        @test g ≈ ifft(fft(g))
    end

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
        @test size(wfk) == (1, 4, 8)
    end

    @testset "VASP outputs" begin
        @test natom(poscar) == 4
        @test atomnames(poscar) == ["Ir", "Si"]
        @test atomtypes(poscar) == [77, 14]
        @test natomtypes(poscar) == 2
    end
end
