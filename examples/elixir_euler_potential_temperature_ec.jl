using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

function initial_condition_density_wave(x, t,
                                        equations::CompressibleEulerPotentialTemperatureEquations1D)
    v1 = 1
    rho = 1 + exp(sinpi(2 * x[1]))
    p = 1
    return prim2cons(SVector(rho, v1, p), equations)
end

equations = CompressibleEulerPotentialTemperatureEquations1D(1004, 717)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)
surface_flux = flux_ec
volume_flux = flux_ec
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0,)
coordinates_max = (1.0,)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
                n_cells_max = 100_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_density_wave, solver)
###############################################################################
# ODE solvers, callbacks etc.
tspan = (0.0, 0.4)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (energy_kinetic,
                                                                 energy_total, entropy,
                                                                 pressure))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

###############################################################################
# run the simulation
sol = solve(ode,
            SSPRK43(thread = Trixi.True());
            maxiters = 1.0e7,
            ode_default_options()..., callback = callbacks)
