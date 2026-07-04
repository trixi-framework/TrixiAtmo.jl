@muladd begin
#! format: noindent

# Calculate 1D flux for a single point
@inline function flux(u, orientation::Integer,
                      equations::CompressibleEulerAtmo{NDIMS}) where {NDIMS}
    normal_direction = SVector{NDIMS}(i == orientation ? 1 : 0 for i in 1:NDIMS)
    return flux(u, normal_direction, equations)
end

@inline function flux(u, normal_direction::AbstractVector,
                      equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS, NCONDENS}) where
                 {NDIMS, NVARS, NGAS, NCONDENS}
    rho_total = density_total(u, equations)
    v_normal = dot(vars_moment(u, equations), normal_direction) / rho_total
    p = pressure(u, equations)

    # momentum equations
    f_mom = vars_moment(u, equations) .* v_normal + normal_direction .* p

    # thermodynamic_equation
    f_td = flux_td(u, equations, equations.td_equation, equations.td_state) * v_normal

    # mass equations
    f_mass_air = vars_airborn(u, equations) .* v_normal
    f_mass_precip1 = vars_precip(u, equations) .* v_normal
    f_mass_passive = vars_passive(u, equations) .* v_normal

    # TODO
    #f_mom_precip, f_td_precip, f_mass_precip2 = flux_precip(u, normal_direction, equations)

    #ret = SVector((f_mom + f_mom_precip)..., f_td + f_td_precip,
    #                f_mass_air..., (f_mass_precip1 + f_mass_precip2)..., f_mass_passive...)

    return SVector(f_mom..., f_td, f_mass_air..., f_mass_precip1..., f_mass_passive...)
end

#=
@inline function flux_precip(u, normal_direction::AbstractVector,
                             equations::CompressibleEulerAtmo{NDIMS, NVARS, NGAS,
                                                              NCONDENS}) where
                 {NDIMS, NVARS, NGAS, NCONDENS}

    # precipitation along third direction
    # TODO this is not correct for curved geometries (earth)
    # TODO earth normal componente instead of normal_direction[NDIMS]

    rho_total = density_total(u, equations)
    v_precip = velocities_precip(u, equations)
    f_mom_z = vars_moment(u, equations)[NDIMS] * normal_direction[NDIMS] / rho_total *
              dot(v_precip, vars_precip(u, equations))
    f_mass_precip = vars_precip(u, equations) .* v_precip * normal_direction[NDIMS]

    return SVector(ntuple(i -> 0, Val(NDIMS - 1))..., f_mom_z),
           0,
           f_mass_precip
end
=#

# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.
@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerAtmo{NDIMS, NVARS}) where {
                                                                                                NDIMS,
                                                                                                NVARS
                                                                                                }
    a = flux_lmars.speed_of_sound
    norm_ = norm(normal_direction)

    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_n_ll = dot(vars_moment(u_ll, equations), normal_direction) / rho_total_ll
    v_n_rr = dot(vars_moment(u_rr, equations), normal_direction) / rho_total_rr
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)

    rho = 0.5f0 * (rho_total_ll + rho_total_rr)

    p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v_n_rr - v_n_ll) / norm_
    v_interface = 0.5f0 * (v_n_ll + v_n_rr) - 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

    if (v_interface >= 0)
        f = u_ll * v_interface
        f_td = flux_lmars_td(p_ll, v_interface, equations, equations.td_equation)
    else
        f = u_rr * v_interface
        f_td = flux_lmars_td(p_rr, v_interface, equations, equations.td_equation)
    end

    # additional terms in momentum equation, pad with zeros
    f_mom = normal_direction * p_interface

    return f + SVector(f_mom...,
                       f_td,
                       ntuple(i -> 0, Val(NVARS - NDIMS - 1))...)
end

"""
    flux_kennedy_gruber(u_ll, u_rr, orientation_or_normal_direction,
                        equations::CompressibleEulerEquations3D)

Kinetic energy preserving two-point flux by
- Kennedy and Gruber (2008)
  Reduced aliasing formulations of the convective terms within the
  Navier-Stokes equations for a compressible fluid
  [DOI: 10.1016/j.jcp.2007.09.020](https://doi.org/10.1016/j.jcp.2007.09.020)
"""
@inline function flux_kennedy_gruber(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleEulerAtmo{NDIMS, NVARS}) where {
                                                                                            NDIMS,
                                                                                            NVARS
                                                                                            }
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_ll = vars_moment(u_ll, equations) ./ rho_total_ll
    v_rr = vars_moment(u_rr, equations) ./ rho_total_rr
    v_dot_n_ll = dot(v_ll, normal_direction)
    v_dot_n_rr = dot(v_rr, normal_direction)
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)
    td_ll = var_td(u_ll, equations) / rho_total_ll
    td_rr = var_td(u_rr, equations) / rho_total_rr

    # Average each factor of products in flux
    rho_total_avg = 0.5f0 * (rho_total_ll + rho_total_rr)
    v_avg = 0.5f0 .* (v_ll + v_rr)
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    td_avg = 0.5f0 * (td_ll + td_rr)

    # Calculate fluxes depending on normal_direction
    f_mass_air = v_dot_n_avg * 0.5f0 .* (vars_airborn(u_ll, equations) +
                  vars_airborn(u_rr, equations))
    f_mass_precip1 = v_dot_n_avg * 0.5f0 .* (vars_precip(u_ll, equations) +
                      vars_precip(u_rr, equations))
    f_mass_passive = v_dot_n_avg * 0.5f0 .* (vars_passive(u_ll, equations) +
                      vars_passive(u_rr, equations))
    f_mass_total = v_dot_n_avg * rho_total_avg

    f_mom = f_mass_total * v_avg + normal_direction * p_avg

    f_td = flux_kennedy_gruber_td(f_mass_total, td_avg, p_avg, v_dot_n_avg,
                                  equations, equations.td_equation)
    return SVector(f_mom...,
                   f_td,
                   f_mass_air..., f_mass_precip1..., f_mass_passive...)
end

"""
    flux_ranocha(u_ll, u_rr, orientation_or_normal_direction,
                 equations::CompressibleEulerEquations2D)

Entropy conserving and kinetic energy preserving two-point flux by
- Hendrik Ranocha (2018)
  Generalised Summation-by-Parts Operators and Entropy Stability of Numerical Methods
  for Hyperbolic Balance Laws
  [PhD thesis, TU Braunschweig](https://cuvillier.de/en/shop/publications/7743)
See also
- Hendrik Ranocha (2020)
  Entropy Conserving and Kinetic Energy Preserving Numerical Methods for
  the Euler Equations Using Summation-by-Parts Operators
  [Proceedings of ICOSAHOM 2018](https://doi.org/10.1007/978-3-030-39647-3_42)
"""
@inline function flux_ranocha(u_ll, u_rr, normal_direction::AbstractVector,
                              equations::CompressibleEulerAtmo{NDIMS, NVARS}) where {
                                                                                     NDIMS,
                                                                                     NVARS
                                                                                     }
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_ll = vars_moment(u_ll, equations) ./ rho_total_ll
    v_rr = vars_moment(u_rr, equations) ./ rho_total_rr
    v_dot_n_ll = dot(v_ll, normal_direction)
    v_dot_n_rr = dot(v_rr, normal_direction)
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)

    gamma_ll = gamma_total(vars_gas(u_ll, equations), vars_condens(u_ll, equations),
                           equations.td_state)
    gamma_rr = gamma_total(vars_gas(u_rr, equations), vars_condens(u_rr, equations),
                           equations.td_state)

    # Compute the necessary mean values
    rho_total_mean = ln_mean(rho_total_ll, rho_total_rr)
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_total_ll * p_rr, rho_total_rr * p_ll)
    v_avg = 0.5f0 .* (v_ll + v_rr)
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    velocity_square_avg = 0.5f0 * dot(v_ll, v_rr)
    inv_gamma_minus_one = 1 / (0.5f0 * (gamma_ll + gamma_rr) - 1)

    # Calculate fluxes depending on normal_direction
    f_mass_air = v_dot_n_avg *
                 ln_mean.(vars_airborn(u_ll, equations),
                          vars_airborn(u_rr, equations))
    f_mass_precip1 = v_dot_n_avg *
                     ln_mean.(vars_precip(u_ll, equations),
                              vars_precip(u_rr, equations))
    f_mass_passive = v_dot_n_avg *
                     ln_mean.(vars_passive(u_ll, equations),
                              vars_passive(u_rr, equations))
    f_mass_total = v_dot_n_avg * rho_total_mean

    f_mom = f_mass_total * v_avg + normal_direction * p_avg

    # TODO: only valid for total energy
    f_td = (f_mass_total * (velocity_square_avg + inv_rho_p_mean * inv_gamma_minus_one)
            +
            0.5f0 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll))

    return SVector(f_mom...,
                   f_td,
                   f_mass_air..., f_mass_precip1..., f_mass_passive...)
end

"""
	flux_ec(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquationsPotentialTemperature2D)

Entropy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_ec(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerAtmo)
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_ll = vars_moment(u_ll, equations) ./ rho_total_ll
    v_rr = vars_moment(u_rr, equations) ./ rho_total_rr
    v_dot_n_ll = dot(v_ll, normal_direction)
    v_dot_n_rr = dot(v_rr, normal_direction)
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)

    rho_theta_ll = var_td(u_ll, equations)
    rho_theta_rr = var_td(u_rr, equations)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_total_ll, rho_total_rr)
    v_avg = 0.5f0 .* (v_ll + v_rr)
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f_mass_air = ln_mean.(vars_airborn(u_ll, equations),
                          vars_airborn(u_rr, equations)) * v_dot_n_avg
    f_mass_precip1 = ln_mean.(vars_precip(u_ll, equations),
                              vars_precip(u_rr, equations)) * v_dot_n_avg
    f_mass_passive = ln_mean.(vars_passive(u_ll, equations),
                              vars_passive(u_rr, equations)) * v_dot_n_avg
    f_mass_total = rho_mean * v_dot_n_avg

    f_mom = f_mass_total * v_avg + p_avg * normal_direction

    f_td = flux_ec_td(f_mass_total, rho_theta_ll, rho_theta_rr, rho_total_ll,
                      rho_total_rr,
                      equations, equations.td_equation)

    return SVector(f_mom...,
                   f_td,
                   f_mass_air..., f_mass_precip1..., f_mass_passive...)
end

"""
	flux_tec(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquationsPotentialTemperature2D)

Total energy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_tec(u_ll, u_rr, normal_direction::AbstractVector,
                          equations::CompressibleEulerAtmo)
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_ll = vars_moment(u_ll, equations) ./ rho_total_ll
    v_rr = vars_moment(u_rr, equations) ./ rho_total_rr
    v_dot_n_ll = dot(v_ll, normal_direction)
    v_dot_n_rr = dot(v_rr, normal_direction)
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)

    rho_theta_ll = var_td(u_ll, equations)
    rho_theta_rr = var_td(u_rr, equations)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_total_ll, rho_total_rr)
    v_avg = 0.5f0 .* (v_ll + v_rr)
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f_mass_air = ln_mean.(vars_airborn(u_ll, equations),
                          vars_airborn(u_rr, equations)) * v_dot_n_avg
    f_mass_precip1 = v_dot_n_avg *
                     ln_mean.(vars_precip(u_ll, equations),
                              vars_precip(u_rr, equations))
    f_mass_passive = v_dot_n_avg *
                     ln_mean.(vars_passive(u_ll, equations),
                              vars_passive(u_rr, equations))
    f_mass_total = v_dot_n_avg * rho_mean

    f_mom = f_mass_total * v_avg + p_avg * normal_direction

    f_td = flux_tec_td(rho_theta_ll, rho_theta_rr, v_dot_n_avg,
                       equations, equations.td_equation)

    return SVector(f_mom...,
                   f_td,
                   f_mass_air..., f_mass_precip1..., f_mass_passive...)
end

"""
	flux_etec(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquationsPotentialTemperature2D)

Entropy and total energy conservative two-point flux by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_etec(u_ll, u_rr, normal_direction::AbstractVector,
                           equations::CompressibleEulerAtmo)
    rho_total_ll = density_total(u_ll, equations)
    rho_total_rr = density_total(u_rr, equations)
    v_ll = vars_moment(u_ll, equations) ./ rho_total_ll
    v_rr = vars_moment(u_rr, equations) ./ rho_total_rr
    v_dot_n_ll = dot(v_ll, normal_direction)
    v_dot_n_rr = dot(v_rr, normal_direction)
    p_ll = pressure(u_ll, equations)
    p_rr = pressure(u_rr, equations)

    rho_theta_ll = var_td(u_ll, equations)
    rho_theta_rr = var_td(u_rr, equations)

    # Compute the necessary mean values
    v_avg = 0.5f0 .* (v_ll + v_rr)
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f_td = flux_etec_td(rho_theta_ll, rho_theta_rr, v_dot_n_avg,
                        equations, equations.td_equation)

    f_mass_air = f_td *
                 ln_mean.(vars_airborn(u_ll, equations) ./ rho_theta_ll,
                          vars_airborn(u_rr, equations) ./ rho_theta_rr)
    f_mass_precip1 = f_td *
                     ln_mean.(vars_precip(u_ll, equations) ./ rho_theta_ll,
                              vars_precip(u_rr, equations) ./ rho_theta_rr)
    f_mass_passive = f_td *
                     ln_mean.(vars_passive(u_ll, equations) ./ rho_theta_ll,
                              vars_passive(u_rr, equations) ./ rho_theta_rr)

    f_mass_total = f_td *
                   ln_mean(rho_total_ll / rho_theta_ll, rho_total_rr / rho_theta_rr)

    f_mom = f_mass_total * v_avg + p_avg * normal_direction

    return SVector(f_mom...,
                   f_td,
                   f_mass_air..., f_mass_precip1..., f_mass_passive...)
end
end # @muladd
