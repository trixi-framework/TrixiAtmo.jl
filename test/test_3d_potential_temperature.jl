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

@testitem "Euler potential temperature 3D: elixir_potential_temperature_baroclinic_instability Souza" setup=[
    Setup,
    EulerPotentialTemperature3D
] tags=[:euler_potential_temperature_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        l2=[
                            0.0005402393238279539,
                            0.01550420011077671,
                            0.01381076289049231,
                            0.012941056324550382,
                            0.1399799820386323
                        ],
                        linf=[
                            0.0046516384374393205,
                            1.1039753387835978,
                            0.6663546823502753,
                            1.127630481159228,
                            1.0013116465979692
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        # Avoid intermediate time stops of the `SaveSolutionCallback`, which
                        # make the adaptive time step sensitive to round-off differences
                        dt=100000,
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
                            0.0005196271201210046,
                            0.014181825905456531,
                            0.01231325412075897,
                            0.012925495519842033,
                            0.13992540482870974
                        ],
                        linf=[
                            0.004229567471452977,
                            1.121979078608027,
                            0.6394800112775647,
                            1.1308637692864063,
                            1.001518609998982
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        # Avoid intermediate time stops of the `SaveSolutionCallback`, which
                        # make the adaptive time step sensitive to round-off differences
                        dt=100000,
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
                            0.0005237192322911491,
                            0.014655087121427119,
                            0.012853488121869804,
                            0.012929752307425608,
                            0.13994217749435772
                        ],
                        linf=[
                            0.004240678176210988,
                            1.1141184103385866,
                            0.6517477337297323,
                            1.1292928534049589,
                            1.0016135249983336
                        ],
                        rtol=8e-8,
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        # Avoid intermediate time stops of the `SaveSolutionCallback`, which
                        # make the adaptive time step sensitive to round-off differences
                        dt=100000,
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
                            0.00013512205511938857,
                            0.00013512205496536123,
                            0.00014321100200751858,
                            0.7578315675919848
                        ],
                        linf=[
                            0.023356470865802503,
                            0.0015824285058342448,
                            0.0015824285029173886,
                            0.0005066641719408516,
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
