@testsnippet KernelAbstractionsExamples begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air", "global_circulation")
end

@testitem "KernelAbstractions backend preference" setup=[Setup] tags=[:kernelabstractions] begin
    @test TrixiAtmo.Trixi._PREFERENCE_THREADING == :kernelabstractions
end

@testitem "KernelAbstractions CPU 3D: elixir_potential_temperature_baroclinic_instability native" setup=[
    Setup,
    KernelAbstractionsExamples
] tags=[:kernelabstractions] begin
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
                        dt=20)
    # The kernels of KernelAbstractions.jl are used on the CPU, also for the stages of the
    # Rosenbrock-W method and the vertical Jacobian
    @test Trixi.trixi_backend(ode.u0) isa TrixiAtmo.KernelAbstractions.CPU
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test real(semi.solver) == Float64
    @test real(semi.solver.basis) == Float64
    @test real(semi.solver.mortar) == Float64
    # TODO: `mesh` is currently not `adapt`ed correctly
    @test real(semi.mesh) == Float64
    @test typeof(semi.equations.c_p) == Float64

    @test ode.u0 isa Array
    @test semi.solver.basis.derivative_matrix isa Array

    @test Trixi.storage_type(semi.cache.elements) === Array
    @test Trixi.storage_type(semi.cache.interfaces) === Array
    @test Trixi.storage_type(semi.cache.boundaries) === Array
    @test Trixi.storage_type(semi.cache.mortars) === Array

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)

    # Check that the total mass and the total mass-weighted potential temperature are
    # conserved up to round-off errors also with the kernels of KernelAbstractions.jl
    for v in (1, 5)
        integral_initial = Trixi.integrate((u, equations) -> u[v], ode.u0, semi)
        integral_final = Trixi.integrate((u, equations) -> u[v], sol.u[end], semi)
        @test isapprox(integral_final, integral_initial, rtol = 1e-13)
    end
end

@testitem "KernelAbstractions CPU 3D: elixir_potential_temperature_baroclinic_instability_turbulence native" setup=[
    Setup,
    KernelAbstractionsExamples
] tags=[:kernelabstractions] begin
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
                        dt=20)
    @test Trixi.trixi_backend(ode.u0) isa TrixiAtmo.KernelAbstractions.CPU
    @test semi.cache.upward_unit_vectors isa Array

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_parabolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end

@testitem "KernelAbstractions CPU 3D: elixir_potential_temperature_baroclinic_instability_boundary_layer native" setup=[
    Setup,
    KernelAbstractionsExamples
] tags=[:kernelabstractions] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability_boundary_layer.jl"),
                        # Expected errors are the same as on CPU
                        l2=[
                            0.003031953993863809,
                            0.05144871072439155,
                            0.052389663549330674,
                            0.012601381438528281,
                            0.8067718681549034
                        ],
                        linf=[
                            0.02174801375592339,
                            0.8297212480730628,
                            1.253966949391037,
                            0.9093281464230247,
                            4.833143481634522
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        dt=20)
    @test Trixi.trixi_backend(ode.u0) isa TrixiAtmo.KernelAbstractions.CPU
    @test semi.cache.upward_unit_vectors isa Array
    @test semi.cache.lowest_level isa Array

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_parabolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end

@testitem "KernelAbstractions CPU 3D: elixir_potential_temperature_baroclinic_instability_turbo native" setup=[
    Setup,
    KernelAbstractionsExamples
] tags=[:kernelabstractions] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                        # Expected errors are exactly the same as with CUDA, see
                        # `test_cuda_3d.jl`
                        l2=[
                            0.0005231487476660257,
                            0.016979615333645174,
                            0.01594720348014934,
                            0.01656181338485904,
                            0.14232904565449617
                        ],
                        linf=[
                            0.004219636334529486,
                            1.1709299053252764,
                            0.6920569567036126,
                            1.1056508753293657,
                            0.9997829208047619
                        ],
                        tspan=(0.0, 20.0), trees_per_cube_face=(2, 2))
    # The kernels of KernelAbstractions.jl are used on the CPU
    @test Trixi.trixi_backend(ode.u0) isa TrixiAtmo.KernelAbstractions.CPU
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test real(semi.solver) == Float64
    @test real(semi.solver.basis) == Float64
    @test real(semi.solver.mortar) == Float64
    # TODO: `mesh` is currently not `adapt`ed correctly
    @test real(semi.mesh) == Float64
    @test typeof(semi.equations.c_p) == Float64

    @test ode.u0 isa Array
    @test semi.solver.basis.derivative_matrix isa Array

    @test Trixi.storage_type(semi.cache.elements) === Array
    @test Trixi.storage_type(semi.cache.interfaces) === Array
    @test Trixi.storage_type(semi.cache.boundaries) === Array
    @test Trixi.storage_type(semi.cache.mortars) === Array

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end
