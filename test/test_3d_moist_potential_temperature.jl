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
                            5.39076369e-04,
                            1.21484925e-02,
                            1.47024703e-02,
                            1.93497442e-03,
                            1.40233395e-01,
                            1.96548625e-04,
                            3.98726407e-06,
                            0.0
                        ],
                        linf=[
                            4.69940250e-03,
                            3.00922823e-01,
                            5.60897855e-01,
                            5.32547808e-02,
                            1.00154859e+00,
                            8.14828279e-04,
                            5.16939786e-05,
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
