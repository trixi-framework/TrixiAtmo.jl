@testsnippet EulerMoistPotentialTemperature3D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "precipitation")
end

@testitem "Euler moist potential temperature 3D: elixir_potential_temperature_kessler_baroclinic_instability" setup=[
    Setup,
    EulerMoistPotentialTemperature3D
] tags=[:euler_moist_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_kessler_baroclinic_instability.jl"),
                        l2=[
                            5.39440906e-04,
                            1.21889662e-02,
                            1.47341136e-02,
                            1.93454086e-03,
                            1.40231681e-01,
                            1.96742753e-04,
                            4.05239402e-06,
                            0.0
                        ],
                        linf=[
                            4.70141023e-03,
                            3.00728241e-01,
                            5.60790553e-01,
                            5.33729582e-02,
                            1.00156281e+00,
                            8.15536419e-04,
                            5.25568040e-05,
                            0.0
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)

    # The dynamics and the microphysics conserve the mass of dry air and, before any rain
    # reaches the lower boundary, the mass of water
    dry_air(u, equations) = u[1] - u[6] - u[7] - u[8]
    water(u, equations) = u[6] + u[7] + u[8]
    @test isapprox(Trixi.integrate(dry_air, sol.u[end], semi),
                   Trixi.integrate(dry_air, ode.u0, semi), rtol = 1e-12)
    @test isapprox(Trixi.integrate(water, sol.u[end], semi),
                   Trixi.integrate(water, ode.u0, semi), rtol = 1e-12)
    # The microphysics has produced cloud water
    u = Trixi.wrap_array(sol.u[end], semi)
    @test maximum(view(u, 7, ..)) > 1e-5
end

@testitem "Euler moist potential temperature 3D: elixir_potential_temperature_kessler_baroclinic_instability without FluxTurbo" setup=[
    Setup,
    EulerMoistPotentialTemperature3D
] tags=[:euler_moist_potential_temperature_3d] begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_kessler_baroclinic_instability.jl"),
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2))
    u_ode_specialized = copy(sol.u[end])

    # Same fluxes, but without the SIMD-vectorized flux differencing volume integral
    # selected by `FluxTurbo`
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_kessler_baroclinic_instability.jl"),
                  volume_flux = flux_kennedy_gruber_souza_etal,
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2))
    u_ode = copy(sol.u[end])

    @test u_ode_specialized ≈ u_ode
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
