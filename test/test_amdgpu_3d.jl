@testsnippet AMDGPU3DExamples begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air", "global_circulation")
end

@testitem "AMDGPU 3D: elixir_potential_temperature_baroclinic_instability_turbo Float64 / AMDGPU" setup=[
    Setup,
    AMDGPU3DExamples
] tags=[:AMDGPU] begin
    # Using AMDGPU inside the testitem since otherwise the bindings are hidden by the anonymous modules
    using AMDGPU
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                        # Expected errors similar to reference on CPU
                        l2=[
                            0.0005231487476660257,
                            0.01697962314394029,
                            0.015947211831717847,
                            0.016561823530854555,
                            0.14232904565452106,
                            98.99596252593312
                        ],
                        linf=[
                            0.004219636209495059,
                            1.170930023762024,
                            0.692057012805924,
                            1.105650798544498,
                            0.9997828853549322,
                            333.1178018152714
                        ],
                        tspan=(0.0, 20.0), trees_per_cube_face=(2, 2),
                        storage_type=ROCArray)
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test real(semi.solver) == Float64
    @test real(semi.solver.basis) == Float64
    @test real(semi.solver.mortar) == Float64
    # TODO: `mesh` is currently not `adapt`ed correctly
    @test real(semi.mesh) == Float64
    @test typeof(semi.equations.c_p) == Float64

    @test ode.u0 isa ROCArray
    @test semi.solver.basis.derivative_matrix isa ROCArray

    @test Trixi.storage_type(semi.cache.elements) === ROCArray
    @test Trixi.storage_type(semi.cache.interfaces) === ROCArray
    @test Trixi.storage_type(semi.cache.boundaries) === ROCArray
    @test Trixi.storage_type(semi.cache.mortars) === ROCArray

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end

@testitem "AMDGPU 3D: elixir_potential_temperature_baroclinic_instability_turbo Float32 / AMDGPU" setup=[
    Setup,
    AMDGPU3DExamples
] tags=[:AMDGPU] begin
    # Using AMDGPU inside the testitem since otherwise the bindings are hidden by the anonymous modules
    using AMDGPU
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                        # Expected errors on AMDGPU, which differ from the ones on the CPU
                        # in `Float32` by up to 3e-4 because of the different order of the
                        # operations
                        l2=Float32[0.0005292007,
                                   0.017053418,
                                   0.016058886,
                                   0.01661327,
                                   0.14386404,
                                   101.35988],
                        linf=Float32[0.004470706,
                                     1.1710502,
                                     0.6952989,
                                     1.10515,
                                     1.0601196,
                                     380.0],
                        tspan=(0.0, 20.0), trees_per_cube_face=(2, 2),
                        RealT_for_test_tolerances=Float32,
                        rtol=1.0f-3,
                        real_type=Float32,
                        storage_type=ROCArray)
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test real(semi.solver) == Float32
    @test real(semi.solver.basis) == Float32
    @test real(semi.solver.mortar) == Float32
    # TODO: `mesh` is currently not `adapt`ed correctly
    @test real(semi.mesh) == Float64
    @test typeof(semi.equations.c_p) == Float32

    @test ode.u0 isa ROCArray
    @test semi.solver.basis.derivative_matrix isa ROCArray

    @test Trixi.storage_type(semi.cache.elements) === ROCArray
    @test Trixi.storage_type(semi.cache.interfaces) === ROCArray
    @test Trixi.storage_type(semi.cache.boundaries) === ROCArray
    @test Trixi.storage_type(semi.cache.mortars) === ROCArray

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end
