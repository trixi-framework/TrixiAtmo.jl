@testsnippet EulerPotentialTemperature2D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air")
end

@testitem "Euler potential temperature 2D: elixir_potential_temperature_inertia_gravity_waves" setup=[
    Setup,
    EulerPotentialTemperature2D
] tags=[:euler_potential_temperature_2d, :upstream] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_potential_temperature_inertia_gravity_waves.jl"),
                        l2=[
                            5.368255320277407e-7,
                            0.00013168228515039588,
                            6.980897717202787e-5,
                            3.700138344100214e-5,
                            8.738544195667257e-12
                        ],
                        linf=[
                            3.7341484546704606e-6,
                            0.0005486985487159757,
                            0.00022641115906779307,
                            0.00014530934379308746,
                            4.3655745685100555e-11
                        ], tspan=(0.0, 1800.0), atol=1e-10)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities). `Trixi.solve!` calls the right-hand side inside the
    # timer section "main loop", so the first call outside of it creates new timer
    # sections, which allocates.
    TrixiAtmo.Trixi.rhs_hyperbolic!(similar(sol.u[end]), sol.u[end], semi, sol.t[end])
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)

    # Check that the total mass and the total mass-weighted potential temperature are
    # conserved up to round-off errors. The explicit part of the scheme is conservative, and
    # the implicit stages are conservative as well since the vertical Jacobian uses the same
    # metric terms as the semidiscretization.
    for v in (1, 4)
        integral_initial = Trixi.integrate((u, equations) -> u[v], ode.u0, semi)
        integral_final = Trixi.integrate((u, equations) -> u[v], sol.u[end], semi)
        @test isapprox(integral_final, integral_initial, rtol = 1e-13)
    end

    # The columns of the vertical Jacobian require a mesh that is not periodic in the
    # vertical direction
    mesh_periodic = P4estMesh((4, 2), polydeg = 3, coordinates_min = coordinates_min,
                              coordinates_max = coordinates_max,
                              periodicity = (true, true))
    semi_periodic = SemidiscretizationHyperbolic(mesh_periodic, equations,
                                                 initial_condition, solver,
                                                 boundary_conditions = boundary_condition_periodic)
    @test_throws ArgumentError TrixiAtmo.create_jacobian_cache(VerticalJacobianTwoSided(),
                                                               semi_periodic,
                                                               Trixi.compute_coefficients(0.0,
                                                                                          semi_periodic))
end

@testitem "Euler potential temperature 2D: elixir_potential_temperature_linear_hydrostatic" setup=[
    Setup,
    EulerPotentialTemperature2D
] tags=[:euler_potential_temperature_2d, :upstream] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "mountain_flow",
                                 "elixir_potential_temperature_linear_hydrostatic.jl"),
                        l2=[
                            2.669332229493225e-6,
                            0.0006518362209690377,
                            0.00011039658133803138,
                            0.00022353351014161769,
                            2.647509013539741e-11
                        ],
                        linf=[
                            2.7590006774769193e-5,
                            0.017497225180171938,
                            0.0018056117333556639,
                            0.001078532265978538,
                            1.7462298274040222e-10
                        ],
                        tspan=(0.0, 360.0), cells_per_dimension=(20, 12),
                        atol=1e-5)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 2D: elixir_potential_temperature_linear_nonhydrostatic" setup=[
    Setup,
    EulerPotentialTemperature2D
] tags=[:euler_potential_temperature_2d, :upstream] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "mountain_flow",
                                 "elixir_potential_temperature_linear_nonhydrostatic.jl"),
                        l2=[
                            4.345825921806206e-7,
                            0.00018887049967917287,
                            0.00012755555800319474,
                            8.576762747892351e-5,
                            2.7230966512302478e-11
                        ],
                        linf=[
                            1.0334839461112466e-5,
                            0.013797647921778733,
                            0.007887457557781355,
                            0.0003272987411264694,
                            1.4551915228366852e-10
                        ], tspan=(0.0, 360.0), cells_per_dimension=(20, 12),
                        atol=1e-5)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 2D: elixir_potential_temperature_schaer_mountain" setup=[
    Setup,
    EulerPotentialTemperature2D
] tags=[:euler_potential_temperature_2d, :upstream] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "mountain_flow",
                                 "elixir_potential_temperature_schaer_mountain.jl"),
                        l2=[
                            0.00010948264917453087,
                            0.09903481134935169,
                            0.09822804596320808,
                            0.004246742748260301,
                            1.8456576434476605e-11
                        ],
                        linf=[
                            0.0017570540164735249,
                            2.3792558995181707,
                            2.354152711825134,
                            0.0764797692902448,
                            8.731149137020111e-11
                        ],
                        tspan=(0.0, 360.0), cells_per_dimension=(20, 12),
                        atol=1e-5)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 2D: elixir_potential_temperature_well_balanced_curvilinear" setup=[
    Setup,
    EulerPotentialTemperature2D
] tags=[:euler_potential_temperature_2d, :upstream] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_well_balanced_curvilinear.jl"),
                        l2=[
                            1.9638272579852512e-8,
                            9.794074039239054e-13,
                            9.90302138162671e-13,
                            5.891481774330958e-6,
                            6.256813355477567e-13
                        ],
                        linf=[
                            5.575406336610911e-8,
                            1.5959531831834447e-11,
                            1.575602141459924e-11,
                            1.6726218973417417e-5,
                            3.637978807091713e-12
                        ], atol=1e-5)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 2D: elixir_potential_temperature_well_balanced_curvilinear with Waruszewski/Souza fluxes" setup=[
    Setup,
    EulerPotentialTemperature2D
] tags=[:euler_potential_temperature_2d, :upstream] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "tests",
                                 "elixir_potential_temperature_well_balanced_curvilinear.jl"),
                        l2=[
                            1.9676517857526736e-8,
                            1.5583712851610544e-7,
                            1.9956396607024943e-7,
                            5.902955362653274e-6,
                            6.256813355477567e-13
                        ],
                        linf=[
                            5.695755511681e-8,
                            9.548396023850076e-7,
                            1.3602645625916122e-6,
                            1.708726631477475e-5,
                            3.637978807091713e-12
                        ], surface_flux=(flux_ec, flux_nonconservative_waruszewski_etal),
                        volume_flux=(flux_etec, flux_nonconservative_souza_etal), atol=1e-5)

    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 2D: elixir_potential_temperature_robert_bubble" setup=[
    Setup,
    EulerPotentialTemperature2D
] tags=[:euler_potential_temperature_2d, :upstream] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_potential_temperature_robert_bubble.jl"),
                        l2=[
                            4.3406554813061024e-5,
                            0.0016924594668553596,
                            0.002272895766783168,
                            0.0025209103450235496
                        ],
                        linf=[
                            0.00026253980321122583,
                            0.00764409768320462,
                            0.0111236926691191,
                            0.013215426049498546
                        ],
                        tspan=(0.0, 1.0), cells_per_dimension=(8, 8))

    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler potential temperature 2D: elixir_potential_temperature_robert_bubble with flux_ec/flux_etec" setup=[
    Setup,
    EulerPotentialTemperature2D
] tags=[:euler_potential_temperature_2d, :upstream] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "buoyancy",
                                 "elixir_potential_temperature_robert_bubble.jl"),
                        l2=[
                            4.342413055592755e-5,
                            0.0016945803933380045,
                            0.002274914518054041,
                            0.0025178880253885075
                        ],
                        linf=[
                            0.00026378589582520817,
                            0.00781773220859705,
                            0.011480660134785742,
                            0.012652115689945731
                        ],
                        tspan=(0.0, 1.0), cells_per_dimension=(8, 8), surface_flux=flux_ec,
                        volume_flux=flux_etec)

    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
