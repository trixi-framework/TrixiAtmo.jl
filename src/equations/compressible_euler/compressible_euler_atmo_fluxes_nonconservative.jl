@muladd begin
#! format: noindent

"""
    flux_nonconservative_waruszewski_etal(u_ll, u_rr, aux_ll, aux_rr,
                                          normal_direction::AbstractVector,
                                          equations::CompressibleEulerAtmo)

Well-balanced gravity term for an isothermal background state
-  Maciej Waruszewski and Jeremy E. Kozdon and Lucas C. Wilcox and Thomas H. Gibson and Francis X. Giraldo (2022)
   Entropy stable discontinuous {G}alerkin methods for balance laws
   in non-conservative form: Applications to the {E}uler equations with gravity
   [DOI: 10.1016/j.jcp.2022.111507](https://doi.org/10.1016/j.jcp.2022.111507)

The well-balancedness on curvilinear coordinates was proven by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_nonconservative_waruszewski_etal(u_ll, u_rr,
                                                       aux_ll, aux_rr,
                                                       normal_direction::AbstractVector,
                                                       equations::CompressibleEulerAtmo{NDIMS,
                                                                                        NVARS,
                                                                                        NGAS,
                                                                                        NCONDENS,
                                                                                        NPRECIP,
                                                                                        NPASSIVE,
                                                                                        1};) where {
                                                                                                    NDIMS,
                                                                                                    NVARS,
                                                                                                    NGAS,
                                                                                                    NCONDENS,
                                                                                                    NPRECIP,
                                                                                                    NPASSIVE
                                                                                                    }
    rho_ll = density_total(u_ll, equations)
    rho_rr = density_total(u_rr, equations)
    rho_avg = ln_mean(rho_ll, rho_rr)

    phi_ll = aux_ll[1]
    phi_rr = aux_rr[1]
    phi_jump = phi_rr - phi_ll

    f_mom = (rho_avg * phi_jump) * normal_direction

    f_td = flux_nonconservative_waruszewski_etal_td(u_ll, u_rr, normal_direction,
                                                    rho_avg, phi_jump,
                                                    equations, equations.td_equation)

    zero_u = zero(eltype(u_ll))

    return SVector(f_mom...,
                   f_td,
                   ntuple(i -> zero_u, Val(NVARS - NDIMS - 1))...)
end

"""
    flux_nonconservative_artiano_etal(u_ll, u_rr, aux_ll, aux_rr,
                                      normal_direction::AbstractVector,
                                      equations::CompressibleEulerAtmo)

Well-balanced gravity term for a constant potential temperature background state by
-  Marco Artiano, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2025)
   Structure-Preserving High-Order Methods for the Compressible Euler Equations
   in Potential Temperature Formulation for Atmospheric Flows
   (https://arxiv.org/abs/2509.10311)
"""
@inline function flux_nonconservative_artiano_etal(u_ll, u_rr,
                                                   aux_ll, aux_rr,
                                                   normal_direction::AbstractVector,
                                                   equations::CompressibleEulerAtmo{NDIMS,
                                                                                    NVARS,
                                                                                    1,
                                                                                    0,
                                                                                    0,
                                                                                    NPASSIVE,
                                                                                    1};) where {
                                                                                                NDIMS,
                                                                                                NVARS,
                                                                                                NPASSIVE
                                                                                                }
    rho_ll = density_total(u_ll, equations)
    rho_rr = density_total(u_rr, equations)

    # Only implemented for dry air. Mixture could results in gamma_ll and gamma_rr.
    rho_avg = stolarsky_mean(rho_ll, rho_rr, equations.td_state.gamma_gas[1])

    phi_ll = aux_ll[1]
    phi_rr = aux_rr[1]
    phi_jump = phi_rr - phi_ll

    f_mom = (rho_avg * phi_jump) * normal_direction

    f_td = flux_nonconservative_artiano_etal_td(u_ll, u_rr, normal_direction,
                                                rho_avg, phi_jump,
                                                equations, equations.td_equation)

    zero_u = zero(eltype(u_ll))

    return SVector(f_mom...,
                   f_td,
                   ntuple(i -> zero_u, Val(NVARS - NDIMS - 1))...)
end

"""
    flux_nonconservative_souza_etal(u_ll, u_rr, aux_ll, aux_rr,
                                    normal_direction::AbstractVector,
                                    equations::CompressibleEulerAtmo)

-  Souza et al.
   The Flux-Differencing Discontinuous {G}alerkin Method Applied to
   an Idealized Fully Compressible Nonhydrostatic Dry Atmosphere
   [DOI: 10.1029/2022MS003527] (https://doi.org/10.1029/2022MS003527)
"""
@inline function flux_nonconservative_souza_etal(u_ll, u_rr,
                                                 aux_ll, aux_rr,
                                                 normal_direction::AbstractVector,
                                                 equations::CompressibleEulerAtmo{NDIMS,
                                                                                  NVARS,
                                                                                  NGAS,
                                                                                  NCONDENS,
                                                                                  NPRECIP,
                                                                                  NPASSIVE,
                                                                                  1};) where {
                                                                                              NDIMS,
                                                                                              NVARS,
                                                                                              NGAS,
                                                                                              NCONDENS,
                                                                                              NPRECIP,
                                                                                              NPASSIVE
                                                                                              }
    rho_ll = density_total(u_ll, equations)
    rho_rr = density_total(u_rr, equations)
    rho_avg = 0.5f0 * (rho_ll + rho_rr)

    phi_ll = aux_ll[1]
    phi_rr = aux_rr[1]
    phi_jump = phi_rr - phi_ll

    f_mom = (rho_avg * phi_jump) * normal_direction

    f_td = flux_nonconservative_souza_etal_td(u_ll, u_rr, normal_direction,
                                              rho_avg, phi_jump,
                                              equations, equations.td_equation)

    zero_u = zero(eltype(u_ll))

    return SVector(f_mom...,
                   f_td,
                   ntuple(i -> zero_u, Val(NVARS - NDIMS - 1))...)
end
end # @muladd
