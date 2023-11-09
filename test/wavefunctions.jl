@testset "Wavefunctions" begin
    i = PlanewaveIndex(1, 2, 3, 4, -5, 6)
    @test size(wavecar) === (1, 6, 12, 13, 13, 13)
    @test axes(wavecar) === (Base.OneTo(1), Base.OneTo(6), Base.OneTo(12), -6:6, -6:6, -6:6)
    wavecar[1, 2, 3, 4, -5, 6] = 420
    @test wavecar[1, 1, 1, 4, 0, 0] == -4.900393f-5 + 1.3592809f-5im
    @test iszero(wavecar[1, 1, 1, 5, 0, 0])
    @test iszero(wavecar[1, 1, 1, 5, 0, 0])
    @test wavecar[1, 1, 1, -4, 0, 0] == -4.900393f-5 + 1.3592809f-5im
    @test wavecar[1, 6, 12, -1, -1, -1] == -0.070487864f0 - 0.1318744f0im
    # Test that getindex and setindex! indices are the same
    @test wavecar[1,2,3,4,-5,6] == 420
    # Periodic indexing
    @test wavecar[1, 2, 3, 13, -13, 0] === wavecar[1, 2, 3, 0, 0, 0]
    @test_throws BoundsError wavecar[1, 7, 12, 0, 0, 0]
    @test wavecar[:,:,:] == wavecar[1, 1:6, 1:12] == wavecar
    @test wavecar[i] == wavecar[CartesianIndex(i)]
    @test wavecar[1, 2, 3] isa ReciprocalDataGrid{3}
    @test vec(wavecar[1, 2, 3].data) == vec(wavecar.data[:, 3, 2, 1])
    # Energies and occupancies
    @test wavecar.energies == energies(EnergyOccupancy(wavecar))
    @test wavecar.occupancies == occupancies(EnergyOccupancy(wavecar))
    @test min_energy(wavecar) == minimum(wavecar.energies)
    @test max_energy(wavecar) == maximum(wavecar.energies)
    @test min_energy(wavecar) == min_energy(EnergyOccupancy(wavecar))
    @test max_energy(wavecar) == max_energy(EnergyOccupancy(wavecar))
    @test min_occupancy(wavecar) == minimum(wavecar.occupancies)
    @test max_occupancy(wavecar) == maximum(wavecar.occupancies)
    @test min_occupancy(wavecar) == min_occupancy(EnergyOccupancy(wavecar))
    @test max_occupancy(wavecar) == max_occupancy(EnergyOccupancy(wavecar))
    @test min_energy(wavecar) <= fermi(wavecar) <= max_energy(wavecar)
end
