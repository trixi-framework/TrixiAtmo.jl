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
#P4estMesh 
function mapping_HOHQ(xi, eta) 
    # transform the mesh from 409.6 m  x 30.0 m to 409600.0 x 30000

    x, y = 1000.0 * xi, 1000.0 * eta 

    return SVector(x, y)
end 


polydeg = 3

mesh_file = joinpath("src/meshes", "agnesi_hill.inp")
mesh_P4est = P4estMesh{2}(mesh_file, polydeg = polydeg, mapping = mapping_HOHQ)

###############################################################################
#
function mapping(xi_, eta_)
    xi = xi_ * 204800.0 + 204800.0  # L = 409600.0 m 
    eta = eta_ * 15000.0 + 15000.0  # H = 30000.0 m

    H = 30000.0 #upper boundary
    a = 16000.0 #half width mountain
    h = 0.016   #height mountain

    topo = (h * a^2)/(a^2 + xi^2)

    x = xi
    y = H * (eta - topo) / (H - topo)
    return SVector(x, y)
end

cells_per_dimension = (32,75) 
# in the paper are 3 different resolutions:
# H1: (128, 300), dx = 3200 m,  dz = 100 m
# H2: (64, 150),  dx = 64000 m, dz = 200 m 
# H3: (32, 75),   dx = 12800 m, dz = 400 m 

mesh_Structured = StructuredMesh(cells_per_dimension, mapping,
                      periodicity = true)

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerEquations2D(1004.0 / 717.0)

initial_condition = initial_condition_agnesi_hill
source_term = source

boundary_condition_HOHQ = Dict(:left => BoundaryConditionDirichlet(initial_condition_agnesi_hill),
                    :right => BoundaryConditionDirichlet(initial_condition_agnesi_hill),
                    :bottom => boundary_condition_slip_wall,
                    :top => boundary_condition_slip_wall)

boundary_condition_Structured = (x_neg = boundary_condition_periodic, 
                    x_pos = boundary_condition_periodic, 
                    y_neg = boundary_condition_slip_wall,
                    y_pos = boundary_condition_slip_wall)

basis = LobattoLegendreBasis(polydeg)

surface_flux = FluxLMARS(340.0)

volume_flux = flux_kennedy_gruber
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

semi = SemidiscretizationHyperbolic(mesh_Structured, equations, initial_condition, solver,
                                    source_terms = source_term,
                                    boundary_conditions = boundary_condition_Structured)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 240000.0)

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

#stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)
                        #stepsize_callback)

###############################################################################
# run the simulation

# H1: dt = 100.0, H2: dt = 200.0, H3: dt = 400.0
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            maxiters = 1.0e7,
            dt = 400.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);

summary_callback()