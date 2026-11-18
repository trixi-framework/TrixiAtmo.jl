# The moist baroclinic wave of DCMIP2016 (test case 161) with the Kessler warm-rain
# microphysics. For optimal results consider increasing the resolution, e.g., to
# `trees_per_cube_face = (16, 8)`.
#
# References:
# - Paul A. Ullrich, Thomas Melvin, Christiane Jablonowski, Andrew Staniforth (2013)
#   A proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores
#   https://doi.org/10.1002/qj.2241
# - Paul A. Ullrich et al. (2016)
#   Dynamical Core Model Intercomparison Project (DCMIP2016) test case document
#   https://doi.org/10.5281/zenodo.1298671

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

equations = CompressibleEulerMoistPotentialTemperatureEquations3D(c_p = 1004.5,
                                                                  c_v = 717.5,
                                                                  gravity = EARTH_GRAVITATIONAL_ACCELERATION,
                                                                  microphysics = KesslerMicrophysics())

initial_condition = initial_condition_moist_baroclinic_instability

# Microphysics and Coriolis force of the rotating Earth
@inline function source_terms_microphysics_coriolis(u, aux, x, t,
                                                    equations::CompressibleEulerMoistPotentialTemperatureEquations3D)
    return (source_terms_microphysics(u, aux, x, t, equations) +
            source_terms_coriolis(u, aux, x, t, equations))
end

boundary_conditions = (; inside = boundary_condition_slip_wall,
                       outside = boundary_condition_slip_wall)

polydeg = 4
surface_flux = flux_lmars_combined
# `FluxTurbo` selects a SIMD-vectorized implementation of the flux differencing volume
# integral, which is equivalent to `volume_flux = flux_kennedy_gruber_souza_etal`.
volume_flux = FluxTurbo(flux_kennedy_gruber_souza_etal)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
trees_per_cube_face = (8, 4)
mesh = P4estMeshCubedSphereTopography(trees_per_cube_face..., EARTH_RADIUS, 30000,
                                      polydeg = polydeg, initial_refinement_level = 0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_microphysics_coriolis,
                                    boundary_conditions = boundary_conditions,
                                    aux_field = aux_field_gravity_earth)

###############################################################################
# ODE solvers, callbacks etc.
T = 10 # 10 days
tspan = (0.0, T * SECONDS_PER_DAY) # time in seconds for 10 days

# Setting `real_type` allows to change the real number type, e.g., to `Float32`.
# This is particularly useful when changing the `storage_type` to a GPU array
# type such as `ROCArray` (AMD) or `CuArray` (NVIDIA CUDA).
# On such backends, `tuning` selects performance settings of the kernels, e.g.,
# `Tuning(flux_differencing_kernel = FullSweep())` for the kernel used for the
# flux differencing volume integral.
semi = trixi_adapt(semi; real_type = nothing, storage_type = nothing,
                   tuning = nothing)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_integrals = ())

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = SECONDS_PER_DAY / 4,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)

###############################################################################
# Use a Runge-Kutta method with automatic (error based) time step size control
# Enable threading of the RK method for better performance on multiple threads
tol = 1e-6
sol = solve(ode,
            SSPRK43(thread = Trixi.Threaded());
            abstol = tol, reltol = tol, ode_default_options()...,
            callback = callbacks)
