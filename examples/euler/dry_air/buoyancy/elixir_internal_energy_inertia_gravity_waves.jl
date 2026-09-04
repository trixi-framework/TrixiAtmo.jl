# This test case is used to compute convergence rates via a linearized solution.
# The setup follows the approach commonly adopted in benchmark studies; therefore,
# a fixed CFL number is employed.
#
# References:
# - Michael Baldauf and Slavko Brdar (2013):
#   "An analytic solution for linear gravity waves in a channel as a test
#   for numerical models using the non-hydrostatic, compressible Euler equations"
#   Q. J. R. Meteorol. Soc., DOI: 10.1002/qj.2105
#   https://doi.org/10.1002/qj.2105
#
# - Maciej Waruszewski, Jeremy E. Kozdon, Lucas C. Wilcox, Thomas H. Gibson,
#   and Francis X. Giraldo (2022):
#   "Entropy stable discontinuous Galerkin methods for balance laws
#   in non-conservative form: Applications to the Euler equations with gravity"
#   JCP, DOI: 10.1016/j.jcp.2022.111507
#   https://doi.org/10.1016/j.jcp.2022.111507
#
# - Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025):
#   "Structure-Preserving High-Order Methods for the Compressible Euler Equations
#   in Potential Temperature Formulation for Atmospheric Flows"
#   https://arxiv.org/abs/2509.10311

using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

equations = CompressibleEulerInternalEnergyEquationsWithGravity2D(c_p = 1004,
                                                                  c_v = 717,
                                                                  gravity = 9.81)

surface_flux = (flux_conservative_es, flux_nonconservative_es)
volume_flux = (flux_conservative_etec, flux_nonconservative_etec)
polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

boundary_conditions = (;
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)
trees_per_dimension = (60, 8)

mesh = P4estMesh(trees_per_dimension, polydeg = polydeg,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 periodicity = (true, false))

initial_condition = initial_condition_gravity_waves
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)
tspan = (0.0, 1800.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (entropy,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        stepsize_callback)

sol = solve(ode,
            SSPRK43(thread = Trixi.Threaded());
            maxiters = 1.0e7,
            dt = 1e-1, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks, adaptive = false)
