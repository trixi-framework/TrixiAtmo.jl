@muladd begin
#! format: noindent

@doc raw"""
    SplitCovariantShallowWaterEquations2D{GlobalCoordinateSystem} <:  
        AbstractCovariantEquations{2, 3, GlobalCoordinateSystem, 3}

Alternative flux formulation of [`CovariantShallowWaterEquations2D`](@ref) given on the 
reference element by 
```math
J \frac{\partial}{\partial t}
\left[\begin{array}{c} h \\ hv^1 \\ hv^2 \end{array}\right] 
+
\frac{\partial}{\partial \xi^1} 
\left[\begin{array}{c} J h v^1 \\ J h v^1 v^1 \\ J h v^2 v^1 \end{array}\right]
+ 
\frac{\partial}{\partial \xi^2} 
\left[\begin{array}{c} J h v^2 \\ J h v^1 v^2 \\ J h v^2 v^2 \end{array}\right] 
+
J \left[\begin{array}{c} 0 \\ \Upsilon^1 \\ \Upsilon^2 \end{array}\right] 
= 
J\left[\begin{array}{c}0 \\ s^1 \\ s^2 \end{array}\right].
```
In the above, the non-conservative differential terms in the momentum equations are given by
```math
\Upsilon^i = \frac{1}{2}hv^j\big(G^{ik}\partial_j v_k - \partial_j v^i\big) 
+ ghG^{ij}\partial_j (h + b),
```
where we allow for a variable bottom topography defined by $h_s$, and the algebraic 
momentum source terms implemented in `source_terms_geometric_coriolis` are given by
```math
s^i = -\frac{1}{2}\big(\Gamma_{jk}^i hv^j v^k - G^{ik}\Gamma_{jk}^lh v^j v_l \big) 
- f JG^{ij}\varepsilon_{jk} hv^k.
```
In the above, we employ the same notation as in [`CovariantShallowWaterEquations2D`](@ref) 
(including summation over repeated indices) and note that the covariant velocity components are given by $v_i = G_{ij} v^j$. To obtain an entropy-conservative scheme with respect to 
the total energy
```math
\eta = \frac{1}{2}h(v_1 v^1 + v_2v^2)  + \frac{1}{2}gh^2 + ghb,
```
this equation type should be used with `volume_flux = (flux_ec, flux_nonconservative_ec)`.
!!! warning "Experimental implementation"
    The use of entropy-stable split-form/flux-differencing formulations for covariant 
    equations is an experimental feature and may change in future releases.
"""
struct SplitCovariantShallowWaterEquations2D{GlobalCoordinateSystem, RealT <: Real} <:
       AbstractCovariantShallowWaterEquations2D{GlobalCoordinateSystem}
    gravity::RealT  # acceleration due to gravity
    rotation_rate::RealT  # rotation rate for Coriolis term 
    global_coordinate_system::GlobalCoordinateSystem
    function SplitCovariantShallowWaterEquations2D(gravity::RealT,
                                                   rotation_rate::RealT;
                                                   global_coordinate_system = GlobalCartesianCoordinates()) where {RealT <:
                                                                                                                   Real}
        return new{typeof(global_coordinate_system), RealT}(gravity, rotation_rate)
    end
end

# This alternative flux formulation has non-conservative terms even in the absence of 
# variable bottom topography
Trixi.have_nonconservative_terms(::SplitCovariantShallowWaterEquations2D) = True()

# Flux as a function of the state vector u, as well as the auxiliary variables aux_vars, 
# which contain the geometric information required for the covariant form.
# Note that this flux does not include the pressure term.
@inline function Trixi.flux(u, aux_vars, orientation::Integer,
                            equations::SplitCovariantShallowWaterEquations2D)
    # Geometric variables
    J = area_element(aux_vars, equations)

    # Physical variables
    h_vcon = momentum_contravariant(u, equations)
    vcon = velocity_contravariant(u, equations)

    # Compute and store mass flux 
    mass_flux = J * h_vcon[orientation]

    # Reuse mass flux to compute momentum flux
    return SVector(mass_flux, mass_flux * vcon[1], mass_flux * vcon[2])
end

@doc raw"""
    Trixi.flux_ec(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                               orientation::Integer,
                               equations::SplitCovariantShallowWaterEquations2D)

Symmetric part of an entropy-conservative flux for the shallow water equations in covariant 
form. Note that this does not include the pressure term or the non-symmetric curvature 
correction term. When used with [`flux_nonconservative_ec`](@ref) for the nonconservative volume and surface terms, this flux recovers the formulation described in the following 
paper for the special case of the Euclidean metric $G_{ab} = \delta_{ab}$:
- N. Wintermeyer, A. R. Winters, G. J. Gassner, and D. A. Kopriva (2017). An entropy stable
  nodal discontinuous Galerkin method for the two dimensional shallow water equations on 
  unstructured curvilinear meshes with discontinuous bathymetry. Journal of Computational 
  Physics 300:240-242. 
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
!!! warning "Experimental implementation"
    The use of entropy-stable split-form/flux-differencing formulations for covariant 
    equations is an experimental feature and may change in future releases.
"""
@inline function Trixi.flux_ec(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                               orientation::Integer,
                               equations::SplitCovariantShallowWaterEquations2D)
    # Geometric variables
    J_ll = area_element(aux_vars_ll, equations)
    J_rr = area_element(aux_vars_rr, equations)

    # Physical variables
    h_vcon_ll = momentum_contravariant(u_ll, equations)
    h_vcon_rr = momentum_contravariant(u_rr, equations)
    vcon_ll = velocity_contravariant(u_ll, equations)
    vcon_rr = velocity_contravariant(u_rr, equations)

    # Mass flux is simple average
    mass_flux = 0.5f0 * (J_ll * h_vcon_ll[orientation] + J_rr * h_vcon_rr[orientation])

    # Momentum flux is average of mass flux times average of velocities
    return SVector(mass_flux, 0.5f0 * (vcon_ll[1] + vcon_rr[1]) * mass_flux,
                   0.5f0 * (vcon_ll[2] + vcon_rr[2]) * mass_flux)
end

@doc raw"""
    flux_nonconservative_ec(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                            orientation::Integer,
                            equations::SplitCovariantShallowWaterEquations2D)

  Non-symmetric part of an entropy-conservative flux for the shallow water equations in 
  covariant form, consisting of pressure and bottom topography terms, as well as a 
  curvature correction term. This can be used for both the volume and surface terms.
!!! warning "Experimental implementation"
    The use of entropy-stable split-form/flux-differencing formulations for covariant 
    equations is an experimental feature and may change in future releases.
"""
@inline function flux_nonconservative_ec(u_ll, u_rr, aux_vars_ll,
                                         aux_vars_rr,
                                         orientation::Integer,
                                         equations::SplitCovariantShallowWaterEquations2D)
    # Geometric variables
    Gcon_ll = metric_contravariant(aux_vars_ll, equations)
    Gcov_rr = metric_covariant(aux_vars_rr, equations)
    J_ll = area_element(aux_vars_ll, equations)
    h_s_jump = bottom_topography(aux_vars_rr, equations) -
               bottom_topography(aux_vars_ll, equations)

    # Physical variables
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    h_vcon_ll = momentum_contravariant(u_ll, equations)
    vcon_rr = velocity_contravariant(u_rr, equations)
    vcov_rr = Gcov_rr * vcon_rr

    geometric_term = 0.5f0 * h_vcon_ll[orientation] * (Gcon_ll * vcov_rr - vcon_rr)
    pressure_term = equations.gravity * Gcon_ll[:, orientation] * h_ll *
                    (h_rr + h_s_jump)

    return SVector(zero(eltype(u_ll)), J_ll * (geometric_term[1] + pressure_term[1]),
                   J_ll * (geometric_term[2] + pressure_term[2]))
end

@doc raw"""
    flux_nonconservative_surface_simplified(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                            orientation::Integer,
                                            equations::SplitCovariantShallowWaterEquations2D)

  For bottom topography which is continuous across element interfaces, we can significantly 
  simplify the nonconservative surface terms, such that only the pressure contribution 
  remains. In such cases, this flux is equivalent to [`flux_nonconservative_ec`](@ref) when used as a surface flux, but should not be used as a volume flux.
!!! warning "Experimental implementation"
    The use of entropy-stable split-form/flux-differencing formulations for covariant 
    equations is an experimental feature and may change in future releases.
"""
@inline function flux_nonconservative_surface_simplified(u_ll, u_rr, aux_vars_ll,
                                                         aux_vars_rr,
                                                         orientation::Integer,
                                                         equations::SplitCovariantShallowWaterEquations2D)
    # Geometric variables
    Gcon_ll = metric_contravariant(aux_vars_ll, equations)
    J_ll = area_element(aux_vars_ll, equations)

    # Physical variables
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)

    pressure_term = equations.gravity * Gcon_ll[:, orientation] * h_ll * h_rr

    return SVector(zero(eltype(u_ll)), J_ll * pressure_term[1], J_ll * pressure_term[2])
end

# Geometric and Coriolis source terms for a rotating sphere for use with the modified 
# split covariant formulation
@inline function source_terms_geometric_coriolis(u, x, t, aux_vars,
                                                 equations::SplitCovariantShallowWaterEquations2D)
    # Geometric variables
    Gcov = metric_covariant(aux_vars, equations)
    Gcon = metric_contravariant(aux_vars, equations)
    (Gamma1, Gamma2) = christoffel_symbols(aux_vars, equations)
    J = area_element(aux_vars, equations)

    # Physical variables
    h_vcon = momentum_contravariant(u, equations)
    v_con = velocity_contravariant(u, equations)

    # Doubly-contravariant and mixed inertial flux tensors
    h_vcon_vcon = h_vcon * v_con'
    h_vcov_vcon = Gcov * h_vcon_vcon

    # Coriolis parameter
    f = 2 * equations.rotation_rate * x[3] / norm(x)  # 2Ωsinθ

    # Geometric source term
    s_geo = 0.5f0 * (SVector(sum(Gamma1 .* h_vcon_vcon), sum(Gamma2 .* h_vcon_vcon)) -
             Gcon * (Gamma1 * h_vcov_vcon[1, :] + Gamma2 * h_vcov_vcon[2, :]))

    # Combined source terms
    source_1 = s_geo[1] + f * J * (Gcon[1, 2] * h_vcon[1] - Gcon[1, 1] * h_vcon[2])
    source_2 = s_geo[2] + f * J * (Gcon[2, 2] * h_vcon[1] - Gcon[2, 1] * h_vcon[2])

    # Do not scale by Jacobian since apply_jacobian! is called before this
    return SVector(zero(eltype(u)), -source_1, -source_2)
end
end # muladd
