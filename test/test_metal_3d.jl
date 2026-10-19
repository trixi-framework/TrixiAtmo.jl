@testsnippet Metal3DExamples begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air", "global_circulation")
end

@testitem "Metal 3D: elixir_potential_temperature_baroclinic_instability_turbo Float32 / Metal" setup=[
    Setup,
    Metal3DExamples
] tags=[:Metal] begin
    # Using Metal inside the testitem since otherwise the bindings are hidden by the anonymous modules
    using Metal
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                        # Expected errors on Metal, which differ from the ones on the CPU
                        # in `Float32` by up to 3e-4 because of the different order of the
                        # operations
                        l2=Float32[0.0005291956,
                                   0.017055795,
                                   0.01606119,
                                   0.016612822,
                                   0.14386341,
                                   101.35988],
                        linf=Float32[0.0044709444,
                                     1.1709689,
                                     0.69529235,
                                     1.1051824,
                                     1.0602112,
                                     380.0],
                        tspan=(0.0, 20.0), trees_per_cube_face=(2, 2),
                        # Metal does not support `Float64`, so the tolerances of the
                        # adaptive time integrator must be `Float32`
                        tol=1.0f-6,
                        RealT_for_test_tolerances=Float32,
                        real_type=Float32,
                        storage_type=MtlArray)
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test real(semi.solver) == Float32
    @test real(semi.solver.basis) == Float32
    @test real(semi.solver.mortar) == Float32
    # TODO: `mesh` is currently not `adapt`ed correctly
    @test real(semi.mesh) == Float64
    @test typeof(semi.equations.c_p) == Float32

    @test ode.u0 isa MtlArray{Float32}
    @test semi.solver.basis.derivative_matrix isa MtlArray

    @test Trixi.storage_type(semi.cache.elements) === MtlArray
    @test Trixi.storage_type(semi.cache.interfaces) === MtlArray
    @test Trixi.storage_type(semi.cache.boundaries) === MtlArray
    @test Trixi.storage_type(semi.cache.mortars) === MtlArray

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end
