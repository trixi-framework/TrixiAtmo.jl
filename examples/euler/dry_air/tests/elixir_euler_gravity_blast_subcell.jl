
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiAtmo

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerEquationsWithGravityNoPressure2D(1.4)

function initial_condition_constant(x, t,
                                    equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho = exp(-x[2])
    v1 = 0.0
    v2 = 0.0
    p = exp(-x[2])
    inicenter = SVector(0, 0.1)
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)

    # Setup based on example 35.1.4 in https://flash.rochester.edu/site/flashcode/user_support/flash4_ug_4p8.pdf
    r0 = 0.21875f0 # = 3.5 * smallest dx (for domain length=4 and max-ref=6)
    # r0 = 0.5 # = more reasonable setup
    E = 1.0
    p0_inner = 3 * (equations.gamma - 1) * E / (3 * pi * r0^2)

    p = r > r0 ? p : p + p0_inner

    prim = SVector(rho, v1, v2, p, x[2])
    return prim2cons(prim, equations)
end

@inline function Trixi.get_boundary_outer_state(u_inner, t,
                                                boundary_condition::typeof(boundary_condition_slip_wall),
                                                orientation, boundary_index,
                                                mesh::TreeMesh{2},
                                                equations::CompressibleEulerEquationsWithGravityNoPressure2D,
                                                dg, cache, indices...)
    if orientation == 1
        #if boundary_index == 1 # Element is on the right, boundary on the left
        return SVector(u_inner[1],
                       u_inner[2] - 2 * u_inner[2],
                       u_inner[3],
                       u_inner[4],
                       u_inner[5])
    else
        return SVector(u_inner[1],
                       u_inner[2],
                       u_inner[3] - 2 * u_inner[3],
                       u_inner[4],
                       u_inner[5])
    end
end

# See Section 2.3 of the reference below for a discussion of robust
# subsonic inflow/outflow boundary conditions.
#
# - Jan-Reneé Carlson (2011)
#   Inflow/Outflow Boundary Conditions with Application to FUN3D.
#   [NASA TM 20110022658](https://ntrs.nasa.gov/citations/20110022658)
@inline function boundary_condition_subsonic(u_inner, orientation::Integer,
                                             direction, x, t,
                                             surface_flux_function,
                                             equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    rho_loc, v1_loc, v2_loc, p_loc, phi_loc = cons2prim(u_inner, equations)

    # For subsonic boundary: Take pressure from initial condition
    p_out = pressure(initial_condition_constant(x, t, equations), equations)

    prim = SVector(rho_loc, v1_loc, v2_loc, p_out, phi_loc)
    u_surface = prim2cons(prim, equations)

    flux_conservative, flux_nonconservative = surface_flux_function

    # get the appropriate normal vector from the orientation
    if orientation == 1
        normal_direction = SVector(one(eltype(u_inner)), zero(eltype(u_inner)))
    else # orientation == 2
        normal_direction = SVector(zero(eltype(u_inner)), one(eltype(u_inner)))
    end

    if isodd(direction)
        return -flux_conservative.numerical_flux(u_inner, u_surface, -normal_direction,
                                                 equations),
               -flux_nonconservative.numerical_flux(u_inner, u_surface, -normal_direction,
                                                    equations)
    else
        return flux_conservative.numerical_flux(u_inner, u_surface, normal_direction,
                                                equations),
               flux_nonconservative.numerical_flux(u_inner, u_surface, normal_direction,
                                                   equations)
    end
end

@inline function Trixi.get_boundary_outer_state(u_inner, t,
                                                boundary_condition::typeof(boundary_condition_subsonic),
                                                orientation, boundary_index,
                                                mesh::TreeMesh{2},
                                                equations::CompressibleEulerEquationsWithGravityNoPressure2D,
                                                dg, cache, indices...)
    rho_loc, v1_loc, v2_loc, p_loc, phi_loc = cons2prim(u_inner, equations)
    x = Trixi.get_node_coords(cache.elements.node_coordinates, equations, dg, indices...)

    # For subsonic boundary: Take pressure from initial condition
    p_out = pressure(initial_condition_constant(x, t, equations), equations)

    prim = SVector(rho_loc, v1_loc, v2_loc, p_out, phi_loc)
    u_surface = prim2cons(prim, equations)
    return u_surface
end

initial_condition = initial_condition_constant

volume_flux = (flux_kennedy_gruber, flux_nonconservative_chandrashekar_isothermal)
surface_flux = (flux_kennedy_gruber, flux_nonconservative_chandrashekar_isothermal)

surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_kennedy_gruber,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction),
                FluxHydrostaticReconstruction(flux_nonconservative_chandrashekar_isothermal,
                                              hydrostatic_reconstruction))

polydeg = 3
basis = LobattoLegendreBasis(polydeg)
limiter_idp = SubcellLimiterIDP(equations, basis;
                                positivity_variables_cons = ["rho"],
                                positivity_variables_nonlinear = [pressure],
                                local_twosided_variables_cons = ["rho"])
volume_integral = VolumeIntegralSubcellLimiting(limiter_idp;
                                                volume_flux_dg = volume_flux,
                                                volume_flux_fv = surface_flux)
#volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_max = (4.0, 7.7)
coordinates_min = (-4.0, -0.3)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 7,
                n_cells_max = 100_000,
                periodicity = (false, false))

boundary_conditions = (;
                       x_neg = boundary_condition_subsonic,
                       x_pos = boundary_condition_subsonic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_subsonic)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_polydeg = polydeg,
                                     save_analysis = true)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim,
                                     extra_node_variables = (:limiting_coefficient,)) #

stepsize_callback = StepsizeCallback(cfl = 0.25)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

stage_callbacks = (SubcellLimiterIDPCorrection(),)

sol = Trixi.solve(ode, Trixi.SimpleSSPRK33(stage_callbacks = stage_callbacks);
                  dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
                  ode_default_options()...,
                  callback = callbacks);
