using OrdinaryDiffEq
using Trixi
using TrixiAtmo: cons2pot, source_terms_rayleigh_sponge, source_terms_gravity

#Hydrostatic test case from J. Simarro, P. Smolikova, J. Vivoda Paper: 
#"An analytical solution of the stationary fully-compressible linear Euler equations over orography" 

# Initial condition
function initial_condition_agnesi_hill(x, t,
                                           equations::CompressibleEulerEquations2D)
    g = 9.81
    c_p = 1004.0 
    c_v = 717.0 
    R = c_p - c_v 

    p_s = 100_000.0
    T_0 = 280.0 

    delta = g / (R*T_0)

    rho_s = p_s * delta / g 

    p_0 = p_s * exp(-delta * x[2])
    rho_0 = rho_s * exp(-delta * x[2])

    v1 = 8.0
    v2 = 0.0
    return prim2cons(SVector(rho_0, v1, v2, p_0), equations)
end

# Source terms 
@inline function source(u, x, t, equations::CompressibleEulerEquations2D)
    return (source_terms_rayleigh_sponge(u, x, t,
                                                       equations::CompressibleEulerEquations2D) +
            source_terms_gravity(u, equations::CompressibleEulerEquations2D))
end #needs slightly different sponge term for the upper boundary

###############################################################################

mesh_file = joinpath("src/meshes", "agnesi_hill.inp")
mesh = P4estMesh{2}(mesh_file, polydeg = polydeg)

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerEquations2D(1004.0 / 717.0)

initial_condition = initial_condition_agnesi_hill
source_term = source

boundary_condition = Dict(:left => boundary_condition_periodic,
                    :right => boundary_condition_periodic,
                    :bottom => boundary_condition_slip_wall,
                    :top => boundary_condition_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

surface_flux = FluxLMARS(340.0)

volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_term,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 120.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
solution_variables = cons2prim

analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:entropy_conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = solution_variables)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            maxiters = 1.0e7,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);

summary_callback()