@testsnippet EulerPotentialTemperature3D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air")
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_taylor_green_vortex" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_taylor_green_vortex.jl"),
                        l2=[
                            0.007135833158819355,
                            0.05920026666636394,
                            0.058302519710783984,
                            0.09239792853460427,
                            0.0026512482044306874
                        ],
                        linf=[
                            0.025211900914064778,
                            0.1626379859810771,
                            0.1567281830166123,
                            0.2914137341434033,
                            0.010285071822802083
                        ],
                        tspan=(0.0, 1.0))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_taylor_green_vortex with flux_ec" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_taylor_green_vortex.jl"),
                        l2=[
                            0.00713580882298077,
                            0.05920027226189781,
                            0.05830253041920023,
                            0.09239789909513721,
                            0.0026512387869924467
                        ],
                        linf=[
                            0.025211360864895616,
                            0.1626378819726536,
                            0.15673002075139453,
                            0.2914131870829597,
                            0.010284993042893986
                        ],
                        tspan=(0.0, 1.0), surface_flux=flux_ec, volume_flux=flux_ec)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_taylor_green_vortex with flux_tec" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_taylor_green_vortex.jl"),
                        l2=[
                            0.007135944983851745,
                            0.05920027807981362,
                            0.05830252313445877,
                            0.09239791128377808,
                            0.0026512477402571037
                        ],
                        linf=[
                            0.02521600221697473,
                            0.16264074629051783,
                            0.15673138446141865,
                            0.29141401520767773,
                            00.010285068999815183
                        ],
                        tspan=(0.0, 1.0), surface_flux=flux_tec, volume_flux=flux_tec)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_taylor_green_vortex with FluxLMARS/flux_tec" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_taylor_green_vortex.jl"),
                        l2=[
                            0.006704842021819753,
                            0.018574880297682488,
                            0.018606626958508878,
                            0.010451109496991842,
                            0.0032071378388946444
                        ],
                        linf=[
                            0.025690924689206418,
                            0.04766269256457434,
                            0.050155417906555394,
                            0.030904595875059383,
                            0.012502671674568866
                        ],
                        tspan=(0.0, 0.2), surface_flux=FluxLMARS(340.0),
                        volume_flux=flux_tec, atol=5e-6)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0030320086112184466,
                            0.05144411394023726,
                            0.05238807358575555,
                            0.012609443441337428,
                            0.8067575925333563,
                            570.8566086563927
                        ],
                        linf=[
                            0.021751907036954288,
                            0.8297225783677513,
                            1.2589734934841315,
                            0.9103791727770297,
                            4.833051748224818,
                            1704.147344700992
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        dt=20)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities). `Trixi.solve!` calls the right-hand side inside the
    # timer section "main loop", so the first call outside of it creates new timer
    # sections, which allocates.
    TrixiAtmo.Trixi.rhs_hyperbolic!(similar(sol.u[end]), sol.u[end], semi, sol.t[end])
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)

    # Check that the total mass and the total mass-weighted potential temperature are
    # conserved up to round-off errors. The explicit part of the scheme is conservative, and
    # the implicit stages are conservative as well since the vertical Jacobian is in flux
    # form with the same quadrature weights and metric terms as the semidiscretization.
    for v in (1, 5)
        integral_initial = Trixi.integrate((u, equations) -> u[v], ode.u0, semi)
        integral_final = Trixi.integrate((u, equations) -> u[v], sol.u[end], semi)
        @test isapprox(integral_final, integral_initial, rtol = 1e-13)
    end

    # The kernels of the vertical Jacobian access the arrays without bounds checks, so the
    # explicit bounds checks before launching them must catch arrays of the wrong size
    jacobian_cache = TrixiAtmo.create_cache_jacobian(jacobian, semi, ode.u0)
    backend = jacobian_cache.backend
    (; u_vertical) = jacobian_cache
    @test isnothing(TrixiAtmo.check_axes_columns(jacobian_cache))
    u_vertical_too_small = similar(u_vertical, Base.front(size(u_vertical))...,
                                   size(u_vertical, 4) - 1)
    @test_throws DimensionMismatch TrixiAtmo.ldiv_interior_forward!(backend,
                                                                    jacobian_cache,
                                                                    u_vertical_too_small,
                                                                    20.0)
    @test_throws DimensionMismatch TrixiAtmo.copy_state_to_columns!(backend,
                                                                    u_vertical_too_small,
                                                                    Trixi.wrap_array(ode.u0,
                                                                                     semi),
                                                                    jacobian_cache, semi)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability Souza" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.003030220313693854,
                            0.05335666293346016,
                            0.05427229406574965,
                            0.012822475231476501,
                            0.8068192103031516,
                            570.8566086563927
                        ],
                        linf=[
                            0.019652674426639427,
                            0.8433669261309379,
                            1.17137194041459,
                            0.9095778771181094,
                            4.8231112405632075,
                            1704.147344700992
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        dt=20,
                        volume_flux=(flux_tec, flux_nonconservative_souza_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities). `Trixi.solve!` calls the right-hand side inside the
    # timer section "main loop", so the first call outside of it creates new timer
    # sections, which allocates.
    TrixiAtmo.Trixi.rhs_hyperbolic!(similar(sol.u[end]), sol.u[end], semi, sol.t[end])
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability FluxTurbo(flux_ec, flux_nonconservative_isothermal)" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_baroclinic_instability.jl"),
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2), dt = 20)
    u_ode_specialized = copy(sol.u[end])

    # Standard implementation with tuples of conservative and nonconservative fluxes
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_baroclinic_instability.jl"),
                  volume_flux = (flux_ec, flux_nonconservative_isothermal),
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2), dt = 20)
    u_ode = copy(sol.u[end])

    @test u_ode_specialized ≈ u_ode
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities). `Trixi.solve!` calls the right-hand side inside the
    # timer section "main loop", so the first call outside of it creates new timer
    # sections, which allocates.
    TrixiAtmo.Trixi.rhs_hyperbolic!(similar(sol.u[end]), sol.u[end], semi, sol.t[end])
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability Artiano" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0030269733595988573,
                            0.0524806653891805,
                            0.05340928705922288,
                            0.012749946942035753,
                            0.8068491880533599,
                            570.8566086563927
                        ],
                        linf=[
                            0.02037252533491385,
                            0.8360193905634581,
                            1.2146162508728153,
                            0.9110979159038155,
                            4.832699312170291,
                            1704.147344700992
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        dt=20,
                        volume_flux=(flux_etec, flux_nonconservative_artiano_etal))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities). `Trixi.solve!` calls the right-hand side inside the
    # timer section "main loop", so the first call outside of it creates new timer
    # sections, which allocates.
    TrixiAtmo.Trixi.rhs_hyperbolic!(similar(sol.u[end]), sol.u[end], semi, sol.t[end])
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

    # Standard implementation with tuples of conservative and nonconservative fluxes
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
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

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability_turbo (Float32)" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                        l2=Float32[0.00052923214,
                                   0.017052146,
                                   0.01605738,
                                   0.016616695,
                                   0.14386448,
                                   101.35988],
                        linf=Float32[0.0044720173,
                                     1.1709296,
                                     0.695263,
                                     1.1052238,
                                     1.0602722,
                                     380.0],
                        tspan=(0.0, 20.0), trees_per_cube_face=(2, 2),
                        RealT_for_test_tolerances=Float32,
                        real_type=Float32)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
    # The elixir uses the adapted `semi` for the ODE problem and the time span
    @test semi === ode.p
    @test eltype(ode.tspan) == Float32
    @test sol.t[end] isa Float32
    @test eltype(sol.u[end]) == Float32
    @test real(semi.solver) == Float32
    @test real(semi.solver.basis) == Float32
    @test real(semi.solver.mortar) == Float32
    @test typeof(semi.equations.c_p) == Float32
    # TODO: `mesh` is currently not `adapt`ed correctly
    @test real(semi.mesh) == Float64
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
