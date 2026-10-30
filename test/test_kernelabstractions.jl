@testsnippet KernelAbstractionsExamples begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air", "global_circulation")
end

@testitem "KernelAbstractions backend preference" setup=[Setup] tags=[:kernelabstractions] begin
    @test TrixiAtmo.Trixi._PREFERENCE_THREADING == :kernelabstractions
end

@testitem "KernelAbstractions CPU 3D: elixir_potential_temperature_baroclinic_instability.jl native" setup=[
    Setup,
    KernelAbstractionsExamples
] tags=[:kernelabstractions] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        # Expected errors similar to reference on CPU
                        l2=[
                            0.06049646764469473,
                            0.8991477386737083,
                            0.8989027811649443,
                            0.032243454466505415,
                            16.378452703937285,
                            12088.812342810517
                        ],
                        linf=[
                            0.4028624079698643,
                            11.767657917694287,
                            11.767656315593566,
                            0.591082133466972,
                            103.08761789192533,
                            42964.04658728838
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        polydeg=3, dt=20.0)
    # The kernels of KernelAbstractions.jl are used on the CPU, also for the stages of the
    # Rosenbrock-W method and the vertical Jacobian
    @test Trixi.trixi_backend(ode.u0) isa TrixiAtmo.KernelAbstractions.CPU
    # The elixir adapts `semi` before creating the ODE problem
    @test semi === ode.p
    @test real(semi.solver) == Float64

    @test ode.u0 isa Array
    @test Trixi.storage_type(semi.cache.elements) === Array
    @test Trixi.storage_type(semi.cache.interfaces) === Array
    @test Trixi.storage_type(semi.cache.boundaries) === Array

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

@testitem "KernelAbstractions CPU 3D: elixir_potential_temperature_baroclinic_instability_turbo.jl native" setup=[
    Setup,
    KernelAbstractionsExamples
] tags=[:kernelabstractions] begin
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
                        tspan=(0.0, 20.0), trees_per_cube_face=(2, 2))
    @test Trixi.trixi_backend(ode.u0) isa TrixiAtmo.KernelAbstractions.CPU
    @test semi === ode.p
    @test real(semi.solver) == Float64

    @test ode.u0 isa Array
    @test Trixi.storage_type(semi.cache.elements) === Array
    @test Trixi.storage_type(semi.cache.interfaces) === Array
    @test Trixi.storage_type(semi.cache.boundaries) === Array

    # Ensure that the RHS computation overwrites existing data in `du` correctly.
    u_ode = copy(ode.u0)
    du_ode = similar(u_ode)
    fill!(du_ode, convert(eltype(du_ode), NaN))
    Trixi.rhs_hyperbolic!(du_ode, u_ode, semi, first(ode.tspan))
    @test all(isfinite, du_ode)
end
