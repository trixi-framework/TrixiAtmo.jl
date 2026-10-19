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
