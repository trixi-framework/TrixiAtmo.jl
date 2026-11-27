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
                            0.0005231487494040615,
                            0.016979615333646277,
                            0.01594720348014793,
                            0.016561813384810754,
                            0.14232904756222564
                        ],
                        linf=[
                            0.004219636334625854,
                            1.1709299053250506,
                            0.6920569567040576,
                            1.1056508753294336,
                            0.9997829208306257
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
    @test semi.cache.aux_vars.aux isa CuArray

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
    @test semi.cache.aux_vars.aux isa CuArray{Float32}

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end

@testitem "CUDA 3D: elixir_potential_temperature_baroclinic_instability Float64 / CUDA" setup=[
    Setup,
    CUDA3DExamples
] tags=[:CUDA] begin
    # Using CUDA inside the testitem since otherwise the bindings are hidden by the anonymous modules
    using CUDA
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        # Expected errors are exactly the same as on CPU
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
                        dt=20,
                        storage_type=CuArray)
    # The stages of the Rosenbrock-W method and the vertical Jacobian run on the GPU
    @test Trixi.trixi_backend(ode.u0) isa CUDABackend
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
    @test semi.cache.aux_vars.aux isa CuArray

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end

@testitem "CUDA 3D: elixir_potential_temperature_baroclinic_instability_turbulence Float64 / CUDA" setup=[
    Setup,
    CUDA3DExamples
] tags=[:CUDA] begin
    # Using CUDA inside the testitem since otherwise the bindings are hidden by the anonymous modules
    using CUDA
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability_turbulence.jl"),
                        # Expected errors are the same as on CPU
                        l2=[
                            0.003031777673411901,
                            0.051523875638834524,
                            0.05246443994827504,
                            0.012594810991886378,
                            0.8067706647208186
                        ],
                        linf=[
                            0.021713721259518204,
                            0.8313541828542199,
                            1.2490179094366172,
                            0.9091888600567941,
                            4.833058699565186
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        dt=20,
                        storage_type=CuArray)
    @test Trixi.trixi_backend(ode.u0) isa CUDABackend
    @test semi.cache.upward_unit_vectors isa CuArray

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_parabolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end

@testitem "CUDA 3D: elixir_potential_temperature_baroclinic_instability Float32 / CUDA" setup=[
    Setup,
    CUDA3DExamples
] tags=[:CUDA] begin
    # Using CUDA inside the testitem since otherwise the bindings are hidden by the anonymous modules
    using CUDA
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        # Expected errors similar to reference on CPU
                        l2=Float32[0.0030331002,
                                   0.051523436,
                                   0.05246645,
                                   0.012602088,
                                   0.807088],
                        linf=Float32[0.021747112,
                                     0.8327257,
                                     1.2596467,
                                     0.91195524,
                                     4.8730164],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        dt=20,
                        RealT_for_test_tolerances=Float32,
                        rtol=1.0f-3,
                        real_type=Float32,
                        storage_type=CuArray)
    # The stages of the Rosenbrock-W method and the vertical Jacobian run on the GPU
    @test Trixi.trixi_backend(ode.u0) isa CUDABackend
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test eltype(ode.tspan) == Float32
    @test sol.t[end] isa Float32
    @test eltype(sol.u[end]) == Float32
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
    @test semi.cache.aux_vars.aux isa CuArray{Float32}

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end
