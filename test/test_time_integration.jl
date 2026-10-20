@testsnippet TimeIntegration begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "dry_air")
end

@testitem "Time integration: SSPRW32 without Jacobian equals SimpleSSPRK33" setup=[
    Setup,
    TimeIntegration
] tags=[:time_integration] begin
    # With `jacobian = nothing`, i.e., a zero Jacobian, the Rosenbrock-W method reduces to
    # its explicit part, the SSPRK(3,3) method of Shu and Osher, which is also available as
    # `Trixi.SimpleSSPRK33`. The right-hand side of the elixir does not depend on the time,
    # so both time integrators yield the same solution up to round-off errors.
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                  tspan = (0.0, 10.0), trees_per_cube_face = (2, 2), sol = nothing)

    # Explicit time step restricted by the vertical acoustic waves
    sol_rosenbrock = Trixi.solve(ode, SSPRW32(); dt = 0.5, callback = callbacks)
    sol_ssprk33 = Trixi.solve(ode, Trixi.SimpleSSPRK33(); dt = 0.5,
                              callback = callbacks)

    @test sol_rosenbrock.t[end] == 10.0
    # Compare the variables separately, since the geopotential (which does not change) is
    # much larger than the other variables and would dominate a comparison of the vectors
    u_rosenbrock = Trixi.wrap_array(sol_rosenbrock.u[end], semi)
    u_ssprk33 = Trixi.wrap_array(sol_ssprk33.u[end], semi)
    u_initial = Trixi.wrap_array(ode.u0, semi)
    for v in eachvariable(equations)
        @test view(u_rosenbrock, v, ..) ≈ view(u_ssprk33, v, ..)
    end
    # The solution changed, i.e., the stages were actually computed
    @test !(view(u_rosenbrock, 2, ..) ≈ view(u_initial, 2, ..))

    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities). `Trixi.solve!` calls the right-hand side inside the
    # timer section "main loop", so the first call outside of it creates new timer
    # sections, which allocates.
    TrixiAtmo.Trixi.rhs_hyperbolic!(similar(sol_rosenbrock.u[end]), sol_rosenbrock.u[end],
                                    semi, sol_rosenbrock.t[end])
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol_rosenbrock, 100)
end

@testitem "Time integration: SSPRW32 with time stops" setup=[
    Setup,
    TimeIntegration
] tags=[:time_integration] begin
    # `SaveSolutionCallback(dt = 0.75)` adds time stops, which shorten the time steps of
    # `dt = 0.5` such that the solution is saved at the multiples of `0.75`. As for the
    # `Trixi.SimpleSSPRK33`, the Rosenbrock-W method without Jacobian must respect them;
    # otherwise, both methods use different time steps and their solutions differ.
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                  tspan = (0.0, 3.0), trees_per_cube_face = (2, 2), sol = nothing)

    save_solution = SaveSolutionCallback(dt = 0.75, save_initial_solution = false,
                                         save_final_solution = false)
    sol_rosenbrock = Trixi.solve(ode, SSPRW32(); dt = 0.5,
                                 callback = CallbackSet(save_solution))
    sol_ssprk33 = Trixi.solve(ode, Trixi.SimpleSSPRK33(); dt = 0.5,
                              callback = CallbackSet(save_solution))

    @test sol_rosenbrock.t[end] == 3.0
    u_rosenbrock = Trixi.wrap_array(sol_rosenbrock.u[end], semi)
    u_ssprk33 = Trixi.wrap_array(sol_ssprk33.u[end], semi)
    for v in eachvariable(equations)
        @test view(u_rosenbrock, v, ..) ≈ view(u_ssprk33, v, ..)
    end
end

@testitem "Time integration: SSPRW32 resize!" setup=[Setup, TimeIntegration] tags=[:time_integration] begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "global_circulation",
                           "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
                  tspan = (0.0, 1.0), trees_per_cube_face = (2, 2), sol = nothing)

    # Resizing the integrator, e.g., for AMR, resizes all vectors of the same length as `u`
    integrator = Trixi.init(ode, SSPRW32(); dt = 0.5)
    new_size = length(ode.u0) + 6
    resize!(integrator, new_size)
    @test length(integrator.u) == length(integrator.du) == length(integrator.u_tmp) ==
          new_size
    @test all(stage -> length(stage) == new_size && iszero(stage), integrator.stages)
end
