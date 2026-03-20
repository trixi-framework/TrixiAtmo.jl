struct FluxHydrostaticReconstruction{NumericalFlux, HydrostaticReconstruction}
    numerical_flux::NumericalFlux
    hydrostatic_reconstruction::HydrostaticReconstruction
end

@inline function (numflux::FluxHydrostaticReconstruction)(u_ll, u_rr,
                                                          orientation_or_normal_direction,
                                                          equations::Trixi.AbstractEquations)
    @unpack numerical_flux, hydrostatic_reconstruction = numflux

    # Create the reconstructed left/right solution states in conservative form
    u_ll_star, u_rr_star = hydrostatic_reconstruction(u_ll, u_rr, equations)

    # Use the reconstructed states to compute the numerical surface flux
    return numerical_flux(u_ll_star, u_rr_star, orientation_or_normal_direction,
                          equations)
end

# Hydrostatic reconstruction from the paper:
# Ziming Chen, Yingjuan Zhang, Gang Li, Shouguo Qian (2022)
# "A well-balanced Runge-Kutta discontinuous Galerkin method for the Euler equations in isothermal 
# hydrostatic state under gravitational field"
# [DOI:10.1016/j.camwa.2022.05.025](https://doi.org/10.1016/j.camwa.2022.05.025)
@inline function hydrostatic_reconstruction(u_ll, u_rr, equations::CompressibleEulerEquationsWithGravityNoPressure2D)
    # Unpack left and right states
    rho_ll, rho_v_ll, rho_e_ll, phi_ll = u_ll
    rho_rr, rho_v_rr, rho_e_rr, phi_rr = u_rr

    # Compute equilibrium potential (# TODO: For general case we need to add phi_num - phi_exact))
    psi_eq_ll = psi_ll # phi_ll - phi_exact_ll
    psi_eq_rr = psi_rr # phi_rr - phi_exact_rr

    # Set RT0 to 1 for now.
    RT0 = 1

    # Compute equlibrium state
    rho_eq_ll = exp(psi_eq_ll / (RT0)) # How do I get T0?   # Replace RT0 with R*T0 once available
    rho_eq_rr = exp(psi_eq_rr / (RT0))
    p_eq_ll = RT0*rho_eq_ll
    p_eq_rr = RT0*rho_eq_rr
    rho_v_eq_ll = 0.0
    rho_v_eq_rr = 0.0
    rho_e_eq_ll = p_eq_ll / (gamma - 1) + rho_eq_ll * phi_ll
    rho_e_eq_rr = p_eq_rr / (gamma - 1) + rho_eq_rr * phi_rr

    u_eq_ll = (rho_eq_ll, rho_v_eq_ll, rho_e_eq_ll, phi_ll)
    u_eq_rr = (rho_eq_rr, rho_v_eq_rr, rho_e_eq_rr, phi_rr)

    # Compute residual contribution
    u_res_ll = u_ll - u_eq_ll
    u_res_rr = u_rr - u_eq_rr

    # Reconstruct phi
    phi_star = max(phi_ll, phi_rr)

    # Compute reconstructed equilibrium potential
    psi_eq_ll = psi_ll + phi_ll - phi_star
    psi_eq_rr = psi_rr + phi_rr - phi_star

    # Compute reconstructed equilibrium state
    rho_eq_ll = exp(psi_eq_ll / (RT0)) # How do I get T0?   # Replace RT0 with R*T0 once available
    rho_eq_rr = exp(psi_eq_rr / (RT0))
    p_eq_ll = RT0*rho_eq_ll
    p_eq_rr = RT0*rho_eq_rr
    rho_v_eq_ll = 0.0
    rho_v_eq_rr = 0.0
    rho_e_eq_ll = p_eq_ll / (gamma - 1) + rho_eq_ll * phi_ll
    rho_e_eq_rr = p_eq_rr / (gamma - 1) + rho_eq_rr * phi_rr

    u_eq_ll = SVector(rho_eq_ll, rho_v_eq_ll, rho_e_eq_ll, phi_ll)
    u_eq_rr = SVector(rho_eq_rr, rho_v_eq_rr, rho_e_eq_rr, phi_rr)

    # Compute reconstructed state
    u_star_ll = u_eq_ll + u_res_ll
    u_star_rr = u_eq_rr + u_res_rr

    return u_star_ll, u_star_rr
end