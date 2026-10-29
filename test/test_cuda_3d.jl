@testsnippet CUDA3DExamples begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air", "global_circulation")
end

@testitem "CUDA 3D: elixir_potential_temperature_baroclinic_instability_turbo Float64 / CUDA" setup=[
    Setup,
    CUDA3DExamples
] tags=[:CUDA] begin
    # Using CUDA inside the testitem since otherwise the bindings are hidden by the anonymous modules
    using CUDA
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                        # Expected errors similar to reference on CPU
                        l2=[
                            0.0005231487476660257,
                            0.01697962311431438,
                            0.015947211723936594,
                            0.016561823433691548,
                            0.14232904565449617,
                            98.99596252593312
                        ],
                        linf=[
                            0.004219636209495059,
                            1.1709300104924298,
                            0.6920569865271409,
                            1.1056507865964287,
                            0.9997828853548185,
                            333.1178018152714
                        ],
                        tspan=(0.0, 20.0), trees_per_cube_face=(2, 2),
                        storage_type=CuArray)
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test real(semi.solver) == Float64
    @test real(semi.solver.basis) == Float64
    @test real(semi.solver.mortar) == Float64
    # TODO: `mesh` is currently not `adapt`ed correctly
    @test real(semi.mesh) == Float64
    @test typeof(semi.equations.c_p) == Float64

    @test ode.u0 isa CuArray
    @test semi.solver.basis.derivative_matrix isa CuArray

    @test Trixi.storage_type(semi.cache.elements) === CuArray
    @test Trixi.storage_type(semi.cache.interfaces) === CuArray
    @test Trixi.storage_type(semi.cache.boundaries) === CuArray
    @test Trixi.storage_type(semi.cache.mortars) === CuArray

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end

@testitem "CUDA 3D: elixir_potential_temperature_baroclinic_instability_turbo Float32 / CUDA" setup=[
    Setup,
    CUDA3DExamples
] tags=[:CUDA] begin
    # Using CUDA inside the testitem since otherwise the bindings are hidden by the anonymous modules
    using CUDA
    @test_trixi_include(joinpath(EXAMPLES_DIR,
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
                        rtol=1.0f-3,
                        real_type=Float32,
                        storage_type=CuArray)
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test real(semi.solver) == Float32
    @test real(semi.solver.basis) == Float32
    @test real(semi.solver.mortar) == Float32
    # TODO: `mesh` is currently not `adapt`ed correctly
    @test real(semi.mesh) == Float64
    @test typeof(semi.equations.c_p) == Float32

    @test ode.u0 isa CuArray
    @test semi.solver.basis.derivative_matrix isa CuArray

    @test Trixi.storage_type(semi.cache.elements) === CuArray
    @test Trixi.storage_type(semi.cache.interfaces) === CuArray
    @test Trixi.storage_type(semi.cache.boundaries) === CuArray
    @test Trixi.storage_type(semi.cache.mortars) === CuArray

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end
