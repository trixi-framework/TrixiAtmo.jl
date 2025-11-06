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

"""
	initial_condition_gravity_waves(x, t,
                                        equations::CompressibleEulerEnergyEquationsWithGravity2D)

Test cases for linearized analytical solution by
-  Baldauf, Michael and Brdar, Slavko (2013)
   An analytic solution for linear gravity waves in a channel as a test 
   for numerical models using the non-hydrostatic, compressible {E}uler equations
   [DOI: 10.1002/qj.2105] (https://doi.org/10.1002/qj.2105)
"""
function initial_condition_gravity_waves(x, t,
                                         equations::CompressibleEulerEnergyEquationsWithGravity2D)
    g = equations.g
    c_p = equations.c_p
    c_v = equations.c_v
    # center of perturbation
    x_c = 100_000.0
    a = 5_000
    H = 10_000
    R = c_p - c_v    # gas constant (dry air)
    T0 = 250
    delta = g / (R * T0)
    DeltaT = 0.001
    Tb = DeltaT * sinpi(x[2] / H) * exp(-(x[1] - x_c)^2 / a^2)
    ps = 100_000  # reference pressure
    rhos = ps / (T0 * R)
    rho_b = rhos * (-Tb / T0)
    p = ps * exp(-delta * x[2])
    rho = rhos * exp(-delta * x[2]) + rho_b * exp(-0.5 * delta * x[2])
    v1 = 20
    v2 = 0
    return prim2cons(SVector(rho, v1, v2, p, g * x[2]), equations)
end

equations = CompressibleEulerEnergyEquationsWithGravity2D(c_p = 1004,
                                                          c_v = 717,
                                                          gravity = 9.81)

# We have an isothermal background state with T0 = 250 K. 
# The reference speed of sound can be computed as:
# cs = sqrt(gamma * R * T0)
cs = sqrt(equations.gamma * equations.R * 250)
surface_flux = (FluxLMARS(cs), flux_zero)
volume_flux = (flux_ranocha, flux_nonconservative_waruzewski_etal)
polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)
cells_per_dimension = (60, 8)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, false))
source_terms = nothing
initial_condition = initial_condition_gravity_waves
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
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
            SSPRK43(thread = Trixi.True());
            maxiters = 1.0e7,
            dt = 1e-1, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks, adaptive = false)
