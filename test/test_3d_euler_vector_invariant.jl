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
                            0.0005185551293504875,
                            0.017027947132866136,
                            0.015090141253565458,
                            0.014702882313897592,
                            0.13991185685924915,
                            98.99596252593312
                        ],
                        linf=[
                            0.004215976505289154,
                            0.9447573700011822,
                            0.5192256997767286,
                            0.9433480065440424,
                            1.0011212827678264,
                            333.1178018152714
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2))
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
