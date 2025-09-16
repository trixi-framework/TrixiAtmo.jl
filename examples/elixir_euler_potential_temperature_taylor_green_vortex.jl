using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

function initial_condition_taylor_green_vortex(x, t,
                                               equations::CompressibleEulerPotentialTemperatureEquations3D)
    rho = 1.0
    v1 = sin(x[1]) * cos(x[2]) * cos(x[3])
    v2 = -cos(x[1]) * sin(x[2]) * cos(x[3])
    v3 = 0.0
    p = 10 +
        1.0 / 16.0 *
        ((cos(2 * x[1]) + cos(2 * x[2]))*(cos(2 * x[1]) + 2) - 2)

    return TrixiAtmo.prim2cons(SVector(rho, v1, v2, v3, p), equations)
end

equations = CompressibleEulerPotentialTemperatureEquations3D()
polydeg = 3
basis = LobattoLegendreBasis(polydeg)
surface_flux = flux_etec
volume_flux = flux_etec
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0, 0.0)
coordinates_max = (1.0, 1.0, 1.0) .* (2 * pi)
cells_per_dimension = (8, 8, 8)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                      periodicity = (true, true, true))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_taylor_green_vortex,
                                    solver)
###############################################################################
# ODE solvers, callbacks etc.
tspan = (0.0, 1.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

###############################################################################
# run the simulation
sol = solve(ode,
            SSPRK43(thread = Trixi.True()),
            maxiters = 1.0e7,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks)
