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
                            0.0005204136042404984,
                            0.015991022109756073,
                            0.014837138450798358,
                            0.014428555186414366,
                            0.1407361147152465,
                            98.99596252593312
                        ],
                        linf=[
                            0.004228248228906217,
                            1.153815702972718,
                            0.6939857758618339,
                            1.1165452201407573,
                            0.9994457574504167,
                            333.1178018152714
                        ],
                        tspan=(0.0, 100.0), trees_per_cube_face=(2, 2),
                        storage_type=CuArray)
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test real(ode.p.solver) == Float64
    @test real(ode.p.solver.basis) == Float64
    @test real(ode.p.solver.mortar) == Float64
    # TODO: `mesh` is currently not `adapt`ed correctly
    @test real(ode.p.mesh) == Float64
    @test typeof(ode.p.equations.c_p) == Float64

    @test ode.u0 isa CuArray
    @test ode.p.solver.basis.derivative_matrix isa CuArray

    @test Trixi.storage_type(ode.p.cache.elements) === CuArray
    @test Trixi.storage_type(ode.p.cache.interfaces) === CuArray
    @test Trixi.storage_type(ode.p.cache.boundaries) === CuArray
    @test Trixi.storage_type(ode.p.cache.mortars) === CuArray

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, ode.p, first(ode.tspan))
    @test all(isfinite, du_ode)
end
