using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

###############################################################################
# Run `elixir_potential_temperature_vortex_shedding.jl`, which saves a restart file every
# `restart_interval` time steps

trees_per_cube_face = (16, 8)
adapt_vertical_grid = GalChen()
t_final = 1.0 * SECONDS_PER_DAY
restart_interval = 1000

trixi_include(@__MODULE__,
              joinpath(@__DIR__, "elixir_potential_temperature_vortex_shedding.jl");
              trees_per_cube_face, adapt_vertical_grid, restart_interval,
              tspan = (0.0, t_final))

###############################################################################
# Restart from the last restart file saved before the middle of the simulation and run
# until `t_final` again, which reproduces the solution of the run above

restart_step = restart_interval * div(sol.stats.naccept, 2 * restart_interval)
restart_filename = joinpath("out", "restart_$(lpad(restart_step, 9, '0')).h5")
mesh = load_mesh(restart_filename)

# `load_mesh` returns a plain `P4estMesh`, which does not know the number of layers of the
# cubed sphere or that `partition!` must keep its columns together on one MPI rank
attach_cubed_sphere_metadata!(mesh; layers = trees_per_cube_face[2],
                              keep_columns_together = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_coriolis_small_earth,
                                    boundary_conditions = boundary_conditions)

tspan = (load_time(restart_filename), t_final)
dt = load_dt(restart_filename)
ode = semidiscretize(semi, tspan, restart_filename)

# Do not overwrite the initial snapshot written by the run above
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)
save_solution = SaveSolutionCallback(dt = 5760.0, save_initial_solution = false,
                                     save_final_solution = true)
callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        save_restart)

integrator = init(ode, SSPRK43(thread = Trixi.Threaded()); dt,
                  abstol = tol, reltol = tol, ode_default_options()...,
                  callback = callbacks)

# Load the context of the time integrator with error-based step size control and the
# number of the last time step
load_adaptive_time_integrator!(integrator, restart_filename)
load_timestep!(integrator, restart_filename)

sol = solve!(integrator)
