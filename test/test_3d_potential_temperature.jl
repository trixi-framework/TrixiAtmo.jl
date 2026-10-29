@testsnippet EulerPotentialTemperature3D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air")
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability Souza" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005402395829626791,
                            0.0155041989645114,
                            0.0138107622096068,
                            0.012941059816206165,
                            0.13997998628314068,
                            98.99596252593278
                        ],
                        linf=[
                            0.00465163790895784,
                            1.103975340100896,
                            0.6663546984155281,
                            1.1276305037303587,
                            1.001311664438731,
                            333.1178018152714
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        atol=4e-9)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability isothermal" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005196272070117161,
                            0.0141817969707941,
                            0.012313242513833094,
                            0.012925522828705729,
                            0.1399254026439438,
                            98.99596252593278
                        ],
                        linf=[
                            0.0042295649514556555,
                            1.121978987805587,
                            0.6394799808438503,
                            1.1308638979672447,
                            1.0015178951935582,
                            333.1178018152714
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        volume_flux=(flux_ec, flux_nonconservative_isothermal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability Artiano" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005237192887527135,
                            0.014655078153251792,
                            0.012853482849295788,
                            0.012929758648242088,
                            0.13994217088888805,
                            98.99596252593278
                        ],
                        linf=[
                            0.00424067775992687,
                            1.1141182921289354,
                            0.6517477285197878,
                            1.1292930037247177,
                            1.0016134087282467,
                            333.1178018152714
                        ],
                        rtol=8e-8,
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        volume_flux=(flux_etec, flux_nonconservative_artiano_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability with combined fluxes" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2))
    u_ode_specialized = copy(sol.u[end])

    # Same combined surface flux, but with the fused combined volume flux instead of the
    # SIMD-vectorized flux differencing volume integral selected by `FluxTurbo`.
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                  volume_flux = flux_kennedy_gruber_souza_etal,
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2))
    u_ode_combined = copy(sol.u[end])

    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_baroclinic_instability.jl"),
                  surface_flux = (FluxLMARS(340), flux_zero),
                  volume_flux = (flux_kennedy_gruber, flux_nonconservative_souza_etal),
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2))
    u_ode = copy(sol.u[end])

    @test u_ode_specialized ≈ u_ode
    @test u_ode_combined ≈ u_ode
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_held_suarez" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    if CI_ON_MACOS
        _rtol = 7e-8  # increased error tolerance
    else
        _rtol = sqrt(eps(Float64))  # default
    end

    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_held_suarez.jl"),
                        l2=[
                            0.0031433373482917877,
                            0.000132274039341287,
                            0.00013227403926677734,
                            0.00014259768410822775,
                            0.7578318727895532,
                            569.8247547308886
                        ],
                        linf=[
                            0.023356419582470034,
                            0.0015226273650469607,
                            0.0015226273524388042,
                            0.0005022233032394358,
                            4.88597072706591,
                            1703.946276059638
                        ],
                        rtol=_rtol,
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY),
                        trees_per_cube_face=(2, 2))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_vortex_shedding" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_vortex_shedding.jl"),
                        l2=[
                            0.00011150504131458794,
                            0.036164709804679646,
                            0.03613391367835746,
                            0.047770193565035125,
                            0.03290767203524783,
                            5.137604831206545
                        ],
                        linf=[
                            0.002989575856775284,
                            0.26758732727117207,
                            0.15524562895330518,
                            0.2569864573155951,
                            0.602254907212,
                            180.54091410903135
                        ],
                        rtol=1e-9,
                        tspan=(0.0, 0.0001 * SECONDS_PER_DAY),
                        trees_per_cube_face=(2, 2))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_vortex_shedding with Sleve" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_vortex_shedding.jl"),
                        l2=[
                            0.00011149963508816154,
                            0.03616430242582638,
                            0.03613385848182371,
                            0.04777013134726074,
                            0.03290700386806988,
                            5.1379057814587235
                        ],
                        linf=[
                            0.0029894693413969886,
                            0.26594803885895973,
                            0.15524563612022302,
                            0.256989875810371,
                            0.6022332266465469,
                            180.54091410903135
                        ],
                        rtol=1e-9,
                        tspan=(0.0, 0.0001 * SECONDS_PER_DAY),
                        trees_per_cube_face=(2, 2), adapt_vertical_grid=Sleve(0.7, 0.8))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
