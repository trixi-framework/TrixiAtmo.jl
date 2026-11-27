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
                        tspan=(0.0, 1.0), surface_flux=(flux_ec, flux_zero),
                        volume_flux=(flux_ec, flux_nonconservative_souza_etal))
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
                        tspan=(0.0, 1.0), surface_flux=(flux_tec, flux_zero),
                        volume_flux=(flux_tec, flux_nonconservative_souza_etal))
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
                        tspan=(0.0, 0.2), surface_flux=(FluxLMARS(340.0), flux_zero),
                        volume_flux=(flux_tec, flux_nonconservative_souza_etal), atol=5e-6)
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
                            0.8067575925333563
                        ],
                        linf=[
                            0.021751907036954288,
                            0.8297225783677513,
                            1.2589734934841315,
                            0.9103791727770297,
                            4.833051748224818
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
        integral_initial = Trixi.integrate((u, aux, equations) -> u[v], ode.u0, semi)
        integral_final = Trixi.integrate((u, aux, equations) -> u[v], sol.u[end], semi)
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
                            0.8068192103031516
                        ],
                        linf=[
                            0.019652674426639427,
                            0.8433669261309379,
                            1.17137194041459,
                            0.9095778771181094,
                            4.8231112405632075
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
                            0.8068491880533599
                        ],
                        linf=[
                            0.02037252533491385,
                            0.8360193905634581,
                            1.2146162508728153,
                            0.9110979159038155,
                            4.832699312170291
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
                        l2=Float32[0.00052934315,
                                   0.017086187,
                                   0.016104799,
                                   0.016569622,
                                   0.14389963],
                        linf=Float32[0.004489064,
                                     1.1706848,
                                     0.6953281,
                                     1.1053133,
                                     1.0652466],
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
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_held_suarez.jl"),
                        l2=[
                            0.0031433358598755116,
                            0.00013512204956949112,
                            0.0001351220494460394,
                            0.00014321100200751858,
                            0.7578315675919848
                        ],
                        linf=[
                            0.023356470865802503,
                            0.0015824283966658955,
                            0.0015824283824796912,
                            0.000506664164275603,
                            4.885984294457614
                        ],
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
                            0.00011150506691663358,
                            0.03616471211270547,
                            0.03613391598611015,
                            0.0477701963652348,
                            0.032907679509536124
                        ],
                        linf=[
                            0.0029895756543576457,
                            0.267587332941324,
                            0.15524568730965205,
                            0.25698655865755315,
                            0.6022548472721496
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
                            0.00011149966069463499,
                            0.03616430473403373,
                            0.03613386078981287,
                            0.047770134147715006,
                            0.03290701134352343
                        ],
                        linf=[
                            0.0029894691400704776,
                            0.265948043433294,
                            0.1552456944822348,
                            0.25698997716322003,
                            0.6022331670256449
                        ],
                        rtol=1e-9,
                        tspan=(0.0, 0.0001 * SECONDS_PER_DAY),
                        trees_per_cube_face=(2, 2), adapt_vertical_grid=Sleve(0.7, 0.8))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
