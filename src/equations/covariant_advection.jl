@muladd begin
#! format: noindent

@doc raw"""
    CovariantLinearAdvectionEquation2D <: AbstractCovariantEquations{2, 3, 3}

A variable-coefficient linear advection equation can be defined on a two-dimensional
manifold $S \subset \mathbb{R}^3$ as
```math
\partial_t h + \nabla_S \cdot (h \vec{v}) = 0,
```
where $\nabla_S \cdot$ is the horizontal divergence operator on $S$. We treat this problem 
as a system of equations in which the first variable is the scalar conserved quantity $h$, 
and the second two are the contravariant components $v^1$ and $v^2$ used in the expansion 
with respect to the covariant basis vectors $\vec{a}_1$ and $\vec{a}_2$ as
```math
\vec{v} = v^1 \vec{a}_1 + v^2 \vec{a}_2.
```
where $\vec{a}_1 = \partial \vec{x} / \partial \xi^1$ and 
$\vec{a}_2 = \partial \vec{x} / \partial \xi^2$ are the so-called covariant basis vectors, 
and $\xi^1$ and $\xi^2$ are the local reference space coordinates. The velocity components 
are spatially varying but assumed to be constant in time, so we do not apply any flux or 
dissipation to such variables. The resulting system is then given on the reference element 
as 
```math
\sqrt{G} \frac{\partial}{\partial t}\left[\begin{array}{c} h \\ v^1 \\ v^2 \end{array}\right] 
+
\frac{\partial}{\partial \xi^1} \left[\begin{array}{c} \sqrt{G} h v^1 \\ 0 \\ 0 \end{array}\right]
+ 
\frac{\partial}{\partial \xi^2} \left[\begin{array}{c} \sqrt{G} h v^2 \\ 0 \\ 0 \end{array}\right] 
= 
\left[\begin{array}{c} 0 \\ 0 \\ 0 \end{array}\right],
```
where $G$ is the determinant of the covariant metric tensor expressed as a matrix with 
entries $G_{ij} =  \vec{a}_i \cdot \vec{a}_j$. Note that the variable advection velocity 
components could also be stored as auxiliary variables, similarly to the geometric 
information.
"""
struct CovariantLinearAdvectionEquation2D <: AbstractCovariantEquations{2, 3, 3} end

function Trixi.varnames(::typeof(cons2cons), ::CovariantLinearAdvectionEquation2D)
    return ("scalar", "vcon1", "vcon2")
end

# Convenience function to extract the velocity
function velocity(u, ::CovariantLinearAdvectionEquation2D)
    return SVector(u[2], u[3])
end

# We will define the "entropy variables" here to just be the scalar variable in the first 
# slot, with zeros in the second and third positions
@inline function Trixi.cons2entropy(u, aux_vars,
                                    equations::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u))
    return SVector{3}(u[1], z, z)
end

# We will define the "entropy variables" here to just be the scalar variable in the first 
# slot, with zeros in the second and third positions
@inline function Trixi.cons2entropy(u, aux_vars,
                                    equations::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u))
    return SVector{3}(u[1], z, z)
end

# Numerical flux plus dissipation for abstract covariant equations as a function of the 
# state vector u, as well as the auxiliary variables aux_vars, which contain the geometric 
# information.
@inline function Trixi.flux(u, aux_vars, orientation::Integer,
                            equations::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u))
    return SVector(area_element(aux_vars, equations) * u[1] * u[orientation + 1], z,
                   z)
end

# Local Lax-Friedrichs dissipation which is not applied to the contravariant velocity 
# components, as they should remain unchanged in time
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr, aux_vars_ll,
                                                              aux_vars_rr,
                                                              orientation_or_normal_direction,
                                                              equations::CovariantLinearAdvectionEquation2D)
    z = zero(eltype(u_ll))
    λ = dissipation.max_abs_speed(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                  orientation_or_normal_direction, equations)
    return -0.5f0 * area_element(aux_vars_ll, equations) * λ *
           SVector(u_rr[1] - u_ll[1], z, z)
end

# Maximum wave speed with respect to the a specific orientation
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                                           orientation::Integer,
                                           equations::CovariantLinearAdvectionEquation2D)
    return max(abs(velocity(u_ll, equations)[orientation]),
               abs(velocity(u_rr, equations)[orientation]))
end

# Maximum wave speeds in each direction for CFL calculation
@inline function Trixi.max_abs_speeds(u, aux_vars,
                                      equations::CovariantLinearAdvectionEquation2D)
    return abs.(velocity(u, equations))
end

# Convert contravariant velocity/momentum components to zonal and meridional components
@inline function contravariant2spherical(u, aux_vars,
                                         equations::CovariantLinearAdvectionEquation2D)
    vlon, vlat = basis_covariant(aux_vars, equations) * velocity(u, equations)
    return SVector(u[1], vlon, vlat)
end

# Convert zonal and meridional velocity/momentum components to contravariant components
@inline function spherical2contravariant(u_spherical, aux_vars,
                                         equations::CovariantLinearAdvectionEquation2D)
    vcon1, vcon2 = basis_covariant(aux_vars, equations) \
                   velocity(u_spherical, equations)
    return SVector(u_spherical[1], vcon1, vcon2)
end
end # @muladd
