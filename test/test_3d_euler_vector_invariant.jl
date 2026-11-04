@testsnippet EulerVectorInvariant3D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air")
end

@testitem "Euler vector invariant 3D: elixir_vector_invariant_baroclinic_instability" setup=[
    Setup,
    EulerVectorInvariant3D
] tags=[:euler_vector_invariant_3d] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "global_circulation",
                                 "elixir_vector_invariant_baroclinic_instability.jl"),
                        l2=[
                            0.0005185551299287267,
                            0.01702794705024047,
                            0.015090140614023183,
                            0.014702877472746486,
                            0.13991185540136739,
                            98.99596252593312
                        ],
                        linf=[
                            0.004215974593237748,
                            0.944757292581601,
                            0.5192256773744095,
                            0.9433480536280751,
                            1.0011207426927058,
                            333.1178018152714
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        # Avoid intermediate time stops of the `SaveSolutionCallback`, which
                        # make the adaptive time step sensitive to round-off differences
                        dt=100000)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Euler vector invariant 3D: elixir_vector_invariant_baroclinic_instability without FluxTurbo" setup=[
    Setup,
    EulerVectorInvariant3D
] tags=[:euler_vector_invariant_3d] begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_vector_invariant_baroclinic_instability.jl"),
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2))
    u_ode_specialized = copy(sol.u[end])

    # Same fluxes, but without the SIMD-vectorized flux differencing volume integral
    # selected by `FluxTurbo`
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_vector_invariant_baroclinic_instability.jl"),
                  volume_flux = flux_etec,
                  tspan = (0.0, 100.0), trees_per_cube_face = (2, 2))
    u_ode = copy(sol.u[end])

    @test u_ode_specialized ≈ u_ode
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
