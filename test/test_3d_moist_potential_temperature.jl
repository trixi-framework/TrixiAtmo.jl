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
                            3.02819246e-03,
                            5.38918469e-02,
                            5.42085921e-02,
                            1.06442312e-02,
                            8.08384841e-01,
                            3.50552220e-04,
                            1.87156755e-05,
                            0.0
                        ],
                        linf=[
                            1.95688143e-02,
                            7.24930048e-01,
                            9.51262545e-01,
                            6.28548893e-02,
                            4.83883063e+00,
                            2.05219668e-03,
                            8.70539884e-05,
                            0.0
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)

    # The dynamics and the microphysics conserve the mass of dry air and, before any rain
    # reaches the lower boundary, the mass of water
    dry_air(u, aux, equations) = u[1] - u[6] - u[7] - u[8]
    water(u, aux, equations) = u[6] + u[7] + u[8]
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
    # Rain water below 5 km and none above, such that the rain falls, evaporates, and
    # leaves the domain through the lower boundary
    function initial_condition_rain(x, t, equations)
        u = initial_condition_moist_baroclinic_instability(x, t, equations)
        aux = aux_field_gravity_earth(x, equations)
        rho, v1, v2, v3, p, q_v, q_c, _ = cons2prim(u, aux, equations)
        z = sqrt(sum(abs2, x)) - EARTH_RADIUS
        q_r = 1.0e-3 * max(1 - z / 5000, 0)^2
        rho_dry = rho / (1 + q_v + q_c)
        return prim2cons(SVector(rho_dry * (1 + q_v + q_c + q_r), v1, v2, v3, p,
                                 q_v, q_c, q_r),
                         aux, equations)
    end

    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_kessler_baroclinic_instability.jl"),
                  initial_condition = initial_condition_rain,
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2))
    u_ode_specialized = copy(sol.u[end])

    # Same fluxes, but without the SIMD-vectorized flux differencing volume integral
    # selected by `FluxTurbo`
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_kessler_baroclinic_instability.jl"),
                  initial_condition = initial_condition_rain,
                  volume_flux = flux_kennedy_gruber_souza_etal,
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2))
    u_ode = copy(sol.u[end])

    @test u_ode_specialized ≈ u_ode
    # The mass of dry air is conserved, while the water decreases since the rain leaves the
    # domain through the lower boundary
    dry_air(u, aux, equations) = u[1] - u[6] - u[7] - u[8]
    water(u, aux, equations) = u[6] + u[7] + u[8]
    @test isapprox(Trixi.integrate(dry_air, sol.u[end], semi),
                   Trixi.integrate(dry_air, ode.u0, semi), rtol = 1e-12)
    @test Trixi.integrate(water, sol.u[end], semi) < Trixi.integrate(water, ode.u0, semi)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
