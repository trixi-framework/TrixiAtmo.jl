# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Vertically implicit Jacobian with the two-sided condensation of the nested
# Schur-complement solver. The interior nodes of each element are condensed by a local
# Schur complement onto the vertical momentum (dense LU factors `schur_element`), and the
# element edges form a band matrix with 7 diagonals (`schur_band`), factorized by a
# banded LU.
# The vertical operator is the one of Eq. (34) of the reference below, transformed by the
# diagonal matrix of the norms `n = |Ja^3|` (see `copy_node_to_columns!` in
# `jacobian_3d_gpu.jl`), with the unit contravariant direction `nu = Ja^3 / n`
# (`compute_vertical_frame`) and the metric factor `n / J` (`compute_metric_factor`).
# It conserves mass and mass-weighted potential temperature.

"""
    VerticalJacobianTwoSided(; c = 340)

Approximate Jacobian of the vertical dynamics (vertically propagating acoustic waves and
gravity) of the [`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`](@ref) on a
[`P4estMeshCubedSphereTopography`](@ref) for implicit time integration methods such as
[`SSPRW32`](@ref), including the nested Schur-complement solver with the two-sided
condensation (banded reduced system) of
- Marco Artiano, Andrés Rueda-Ramírez, Oswald Knoth, Hendrik Ranocha, Gregor Gassner (2026)
  GPU-Accelerated Nested Schur-Complement Solvers for Atmospheric Flow Simulations using
  Horizontally-Explicit, Vertically-Implicit Schemes and Discontinuous Galerkin Methods
  [arXiv: 2610.15775](https://arxiv.org/abs/2610.15775)

Pass it as keyword argument `jacobian` to `Trixi.solve` with an implicit time integration
method, e.g., `Trixi.solve(ode, SSPRW32(); dt, jacobian = VerticalJacobianTwoSided())`. Its
arrays are created by `Trixi.init` for the semidiscretization of the ODE problem, with the
array type of the solution, so that it runs on the same backend (CPU or GPU) as the
semidiscretization. The Jacobian decouples into independent linear systems for each
vertical column of nodes. The dissipation of the vertical acoustic waves is based on the
constant speed of sound `c`.

### Permutation
For a column of `nz` elements of polynomial degree `polydeg` (`M = polydeg + 1`
LGL nodes), the vertical unknowns `(rho, rho*theta, rho*w)` are grouped, by
variable, into a small *local* block per element and a small *global* (edge)
block shared by the whole column:

* All `rho` DOFs of element 1 (nodes `1..M`), then
* The `rho*theta` DOFs `2..M-1` of element 1 (interior only), then
* The `rho*w` DOFs `2..M-1` of element 1 (interior only), then
* ... (repeated, element by element) ...
* Finally, 4 edge DOFs per element, ordered element by element:
  `(rho*theta)_1, (rho*w)_1, (rho*w)_M, (rho*theta)_M`.

Density is *never* condensed: it keeps all `M` nodes, including both edges,
in the local block. Only `rho*theta` and `rho*w` are condensed, and only at
the two edge nodes of each element -- that is what "two-sided" refers to, as
opposed to `JacobianCacheOneSidedPointwise`, which condenses all three
variables but only at one edge.

### Sparsity pattern
Example for `nz = 3`, `polydeg = 3` (`M = 4`), in the permuted ordering above:
⠑⢄⠀⣿⠀⠀⠀⠀⠀⠀⠀⠀⢹⣇⣀⠀⠀⠀
⠠⡀⣵⢟⠀⠀⠀⠀⠀⠀⠀⠀⡜⢣⠀⠀⠀⠀
⠀⠀⠀⠀⠑⢄⠀⣿⠀⠀⠀⠀⠀⠉⢹⣇⣀⠀
⠀⠀⠀⠀⠠⡀⣵⢟⠀⠀⠀⠀⠀⠀⡜⢣⠀⠀ = ⎡A B⎤
⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⣿⠀⠀⠀⠉⢹⣇   ⎣C D⎦
⠀⠀⠀⠀⠀⠀⠀⠀⠠⡀⣵⢟⠀⠀⠀⠀⡜⢣
⠂⠠⠶⣉⠀⠀⠀⠀⠀⠀⠀⠀⢟⣵⣤⠀⠀⠀
⠀⠀⠀⠀⠂⠠⠶⣉⠀⠀⠀⠀⠀⠛⢟⣵⣤⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠂⠠⠶⣉⠀⠀⠀⠛⢟⣵

`A` (block-diagonal, one block per element) holds `rho` (diagonal, plus
gravity's off-diagonal coupling to `rho*w`) and the interior `rho*theta`/
`rho*w` (diagonal, densely coupled to each other via the LGL derivative
operator). `D` (bottom-right) is block-tridiagonal before elimination: each
element's 4 edge DOFs couple to the next element's through the single-valued
LMARS flux. `B`/`C` connect local DOFs to both edges of `rho*theta`/`rho*w`
(dense LGL coupling) and, pointwise, `rho*w`'s edges to the local `rho`.

### What `update_jacobian!`/`solve_jacobian!` do
1. **Local Schur complement** (per element, once per step): eliminate the
   diagonal `rho`/`rho*theta` analytically, leaving a dense
   `(M-2)x(M-2)` complement `SA` on interior `rho*w`, LU-factorized once.
2. **Global Schur complement** (per column, once per step): eliminating every
   element's local DOFs against `SA` leaves a `4*nz x 4*nz` banded system
   `SchurBand` (3 sub-/3 super-diagonals) on the edge DOFs, factorized with a
   generic banded LU -- unlike `JacobianCacheOneSidedPointwise`, which solves
   its smaller Schur complement with a specialized block-tridiagonal (block
   Thomas) algorithm instead.
3. **Per-stage solve**: gather the stage right-hand side, eliminate locally
   against `SA`, forward/backward-substitute against the factorized banded
   `SchurBand` for the edge solution, back-substitute into the local DOFs,
   and scatter back into the DGSEM array.

Every kernel is dispatched via KernelAbstractions.jl to a CPU loop or a GPU
kernel. Columns are independent, so the sequential recurrences (fill, Schur
assembly, banded factor/solve) use one thread per column; gather/scatter use
one thread per node; independent per-element work (forward elimination,
back-substitution) uses one thread per `(column, element)`.
"""
struct VerticalJacobianTwoSided{RealT <: Real}
    speed_of_sound::RealT
end

VerticalJacobianTwoSided(; c = 340) = VerticalJacobianTwoSided(c)

# Arrays of the `VerticalJacobianTwoSided` for a semidiscretization, created by
# `create_jacobian_cache`
struct VerticalJacobianTwoSidedCache{NLAYERS, NNODES, RealT <: Real, BackendType,
                                     MatrixRealT <: DenseMatrix{RealT},
                                     ArrayRealT3 <: DenseArray{RealT, 3},
                                     ArrayRealT4 <: DenseArray{RealT, 4},
                                     ArrayRealT5 <: DenseArray{RealT, 5},
                                     VectorInt <: DenseVector{Int},
                                     DerivativeMatrix <: AbstractMatrix{RealT}}
    ncolumns::Int
    speed_of_sound::RealT
    # Inverse quadrature weight `1 / w_1 = 1 / w_M` of the edge nodes; the weights of the
    # volume terms are part of `derivative_hat`
    inv_boundary_weight::RealT
    backend::BackendType # `KernelAbstractions.Backend` of the arrays
    # Unit contravariant direction `nu` of the vertical reference coordinate (slots 1:3)
    # and normalized norm `n = |Ja^3|` (slot 4), laid out `[slot, i, j, k, element]`
    vertical_frame::ArrayRealT5
    # Metric factor `n / J`, `[column, node, layer]`
    metric_factor::ArrayRealT3
    # Pointwise gravity `g (nu . r_hat)`, `[column, node, layer]`
    gravity::ArrayRealT3
    # Column and layer of every element to copy data to and from the columns
    horizontal_index_of_element::VectorInt
    vertical_index_of_element::VectorInt
    derivative_hat::DerivativeMatrix
    # Vertical state with the variables `1 => rho`, `2 => rho_w` (vertical momentum), and
    # `3 => rho_theta`, `[column, v, node, layer]`; also used as right-hand side
    # and solution of the linear system
    u_vertical::ArrayRealT4
    # Potential temperature and `dp/d(rho_theta)`, `[column, node, layer]`
    theta::ArrayRealT3
    dpdrhotheta::ArrayRealT3
    # Entries of the coupling blocks between neighboring elements,
    # `[column, entry, layer]` and `[column, layer]`
    B1m_34::ArrayRealT3
    B1_1::MatrixRealT
    B1_4::MatrixRealT
    B1p_12::ArrayRealT3
    # LU factors of the local Schur complements `S_e` of the elements in the interior
    # nodes, `[column, node, node, layer]`
    schur_element::ArrayRealT4
    # Global Schur complement `S`, a band matrix with 7 diagonals laid out
    # `[column, diagonal, row]`, and the right-hand side and solution of the condensed
    # system, `[column, row]`
    schur_band::ArrayRealT3
    schur_rhs::MatrixRealT
end

@inline function number_of_layers(::VerticalJacobianTwoSidedCache{NLAYERS}) where {NLAYERS}
    return NLAYERS
end

@inline function number_of_nodes(::VerticalJacobianTwoSidedCache{NLAYERS,
                                                                 NNODES}) where {
                                                                                 NLAYERS,
                                                                                 NNODES}
    return NNODES
end

function create_jacobian_cache(jacobian::VerticalJacobianTwoSided, semi, u_ode)
    (; mesh, equations, solver, cache) = semi
    if !(equations isa CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
        throw(ArgumentError("`VerticalJacobianTwoSided` requires the " *
                            "`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`"))
    end
    nlayers = number_of_layers(mesh)
    nnodes_ = nnodes(solver)
    nelements_ = nelements(solver, cache)
    if nelements_ != Trixi.ntrees(mesh)
        throw(ArgumentError("`VerticalJacobianTwoSided` requires a mesh " *
                            "without refinement"))
    end
    ncolumns = nnodes_^2 * nelements_ ÷ nlayers
    RealT = eltype(cache.elements.node_coordinates)

    # The metric terms are computed on the CPU and moved to the storage type of the
    # semidiscretization
    storage = Trixi.storage_type(cache.elements)
    vertical_indices = compute_vertical_mapping(nlayers, solver, cache)
    vertical_line_indices, vertical_element_indices = vertical_indices
    index_maps = compute_element_index_maps(vertical_indices..., solver, cache)
    horizontal_index_of_element = Adapt.adapt(storage, index_maps[1])
    vertical_index_of_element = Adapt.adapt(storage, index_maps[2])
    vertical_frame = Adapt.adapt(storage, compute_vertical_frame(semi))
    metric_factor = compute_metric_factor(ncolumns, nlayers, semi,
                                          vertical_line_indices,
                                          vertical_element_indices)
    metric_factor = Adapt.adapt(storage, metric_factor)
    gravity = compute_gravity(ncolumns, nlayers, semi, vertical_line_indices,
                              vertical_element_indices)
    gravity = Adapt.adapt(storage, gravity)
    inv_boundary_weight = solver.basis.inverse_weights[1]
    derivative_hat = solver.basis.derivative_hat

    # The remaining arrays are allocated with the backend of the semidiscretization
    backend = get_backend(cache.elements.node_coordinates)
    u_vertical = KernelAbstractions.zeros(backend, RealT, ncolumns, 3, nnodes_, nlayers)
    theta = KernelAbstractions.zeros(backend, RealT, ncolumns, nnodes_, nlayers)
    dpdrhotheta = KernelAbstractions.zeros(backend, RealT, ncolumns, nnodes_, nlayers)
    B1m_34 = KernelAbstractions.zeros(backend, RealT, ncolumns, 2, nlayers)
    B1_1 = KernelAbstractions.zeros(backend, RealT, ncolumns, nlayers)
    B1_4 = KernelAbstractions.zeros(backend, RealT, ncolumns, nlayers)
    B1p_12 = KernelAbstractions.zeros(backend, RealT, ncolumns, 2, nlayers)
    schur_element = KernelAbstractions.zeros(backend, RealT, ncolumns, nnodes_ - 2,
                                             nnodes_ - 2,
                                             nlayers)
    schur_band = KernelAbstractions.zeros(backend, RealT, ncolumns, 7, 4 * nlayers)
    schur_rhs = KernelAbstractions.zeros(backend, RealT, ncolumns, 4 * nlayers)

    return VerticalJacobianTwoSidedCache{nlayers, nnodes_, RealT, typeof(backend),
                                         typeof(B1_1), typeof(gravity),
                                         typeof(schur_element),
                                         typeof(vertical_frame),
                                         typeof(horizontal_index_of_element),
                                         typeof(derivative_hat)}(ncolumns,
                                                                 convert(RealT,
                                                                         jacobian.speed_of_sound),
                                                                 inv_boundary_weight,
                                                                 backend,
                                                                 vertical_frame,
                                                                 metric_factor,
                                                                 gravity,
                                                                 horizontal_index_of_element,
                                                                 vertical_index_of_element,
                                                                 derivative_hat,
                                                                 u_vertical,
                                                                 theta, dpdrhotheta,
                                                                 B1m_34,
                                                                 B1_1, B1_4, B1p_12,
                                                                 schur_element,
                                                                 schur_band, schur_rhs)
end

# Assemble and factorize `I / gamma_dt - J` from the state `u_ode`, see `jacobian.jl`
function update_jacobian!(jacobian::VerticalJacobianTwoSidedCache, u_ode, semi,
                          gamma_dt)
    (; backend, u_vertical) = jacobian
    u = wrap_array(u_ode, semi)

    @trixi_timeit_ext backend timer() "copy state" copy_state_to_columns!(backend,
                                                                          u_vertical,
                                                                          u, jacobian,
                                                                          semi)
    @trixi_timeit_ext backend timer() "fill Jacobian" fill_jacobian!(backend, jacobian,
                                                                     semi.equations,
                                                                     gamma_dt)
    @trixi_timeit_ext backend timer() "boundary Schur complement" schur_complement_boundary!(backend,
                                                                                             jacobian,
                                                                                             gamma_dt)
    @trixi_timeit_ext backend timer() "banded LU" lu_band!(backend, jacobian)

    return nothing
end

# Solve `(I / gamma_dt - J) x = b_ode` in place, see `jacobian.jl`
function solve_jacobian!(b_ode, jacobian::VerticalJacobianTwoSidedCache, semi, gamma_dt)
    (; backend, u_vertical) = jacobian
    b = wrap_array(b_ode, semi)

    @trixi_timeit_ext backend timer() "copy right-hand side" copy_rhs_to_columns!(backend,
                                                                                  u_vertical,
                                                                                  b,
                                                                                  jacobian,
                                                                                  semi)

    # Right-hand side of the condensed system, `schur_rhs = b_D - C * A^{-1} * b_A`
    @trixi_timeit_ext backend timer() "forward substitution" ldiv_interior_forward!(backend,
                                                                                    jacobian,
                                                                                    u_vertical,
                                                                                    gamma_dt)

    @trixi_timeit_ext backend timer() "banded solve" ldiv_band!(backend, jacobian)

    # Solution of the interior degrees of freedom from the boundary ones
    @trixi_timeit_ext backend timer() "backward substitution" ldiv_interior_backward!(backend,
                                                                                      jacobian,
                                                                                      u_vertical,
                                                                                      gamma_dt)

    @trixi_timeit_ext backend timer() "copy solution" copy_solution_from_columns!(backend,
                                                                                  b,
                                                                                  u_vertical,
                                                                                  jacobian,
                                                                                  semi,
                                                                                  gamma_dt)

    return nothing
end

###############################################################################
# Fill of the Jacobian: the potential temperature and `dp/d(rho_theta)` per node, the LU
# factors `schur_element` of the local Schur complements, the coupling entries between
# neighboring elements, and the direct contributions to the boundary Schur complement
# `schur_band`

# Fill of element `layer` of the column `column`; the layers of a column are filled in
# order by the same thread
Base.@propagate_inbounds function fill_jacobian_element!(layer, column, nlayers,
                                                         ::Val{NNODES}, theta_element,
                                                         dpdrhotheta_element,
                                                         schur_correction,
                                                         B1m_34, B1_1, B1_4,
                                                         B1p_12, schur_element,
                                                         schur_band,
                                                         u_vertical,
                                                         metric_factor,
                                                         derivative_hat, speed_of_sound,
                                                         inv_speed_of_sound,
                                                         inv_boundary_weight, dp_factor,
                                                         dp_exponent,
                                                         inv_gamma_dt, gamma_dt,
                                                         gravity, theta,
                                                         dpdrhotheta) where {NNODES}
    RealT = eltype(u_vertical)
    zero_ = zero(RealT)
    # Index of the first of the 4 boundary unknowns of this layer in the condensed system
    index_boundary = (layer - 1) * 4 + 1
    for i in 1:NNODES
        theta_element[i] = u_vertical[column, 3, i, layer] /
                           u_vertical[column, 1, i, layer]
        dpdrhotheta_element[i] = dp_factor *
                                 u_vertical[column, 3, i, layer]^dp_exponent
        theta[column, i, layer] = theta_element[i]
        dpdrhotheta[column, i, layer] = dpdrhotheta_element[i]
    end

    # Local Schur complement; the derivative blocks `A13`, `A23`, and `A32` are computed
    # from `derivative_hat`, `theta_element`, `dpdrhotheta_element`, and `metric_factor`
    for i in 1:(NNODES - 2)
        gravity_i = gravity[column, i + 1, layer]
        for j in 1:(NNODES - 2)
            schur_correction[i, j] = zero_
            for k in 1:(NNODES - 2)
                a32 = metric_factor[column, i + 1, layer] *
                      derivative_hat[i + 1, k + 1] *
                      dpdrhotheta_element[k + 1]
                a23 = metric_factor[column, k + 1, layer] *
                      derivative_hat[k + 1, j + 1] *
                      theta_element[j + 1]
                schur_correction[i, j] += a32 * a23
            end
            a13 = metric_factor[column, i + 1, layer] *
                  derivative_hat[i + 1, j + 1]
            if i == j
                schur_element[column, i, j, layer] = inv_gamma_dt -
                                                     gravity_i * gamma_dt * a13 -
                                                     gamma_dt * schur_correction[i, j]
            else
                schur_element[column, i, j, layer] = -gravity_i * gamma_dt * a13 -
                                                     gamma_dt * schur_correction[i, j]
            end
        end
    end
    lu_dense!(layer, column, schur_element, Val(NNODES - 2))

    # `B1_1` of the lowest and `B1_4` of the highest layer are zero
    if layer > 1
        B1m_34[column, 1, layer] = -0.5f0 *
                                   metric_factor[column, 1, layer] *
                                   inv_boundary_weight
        dpdrhotheta_m = dp_factor *
                        u_vertical[column, 3, NNODES, layer - 1]^dp_exponent
        B1m_34[column, 2, layer] = -dpdrhotheta_m * inv_speed_of_sound *
                                   inv_boundary_weight * 0.5f0 *
                                   metric_factor[column, 1, layer]

        B1_1[column, layer] = dpdrhotheta_element[1] * inv_speed_of_sound * 0.5f0 *
                              metric_factor[column, 1, layer] *
                              inv_boundary_weight
    end

    if layer < nlayers
        B1_4[column, layer] = dpdrhotheta_element[NNODES] * inv_speed_of_sound * 0.5f0 *
                              metric_factor[column, NNODES, layer] *
                              inv_boundary_weight

        dpdrhotheta_p = dp_factor *
                        u_vertical[column, 3, 1, layer + 1]^dp_exponent
        B1p_12[column, 1, layer] = -dpdrhotheta_p * inv_speed_of_sound *
                                   inv_boundary_weight * 0.5f0 *
                                   metric_factor[column, NNODES, layer]
        B1p_12[column, 2, layer] = inv_boundary_weight * 0.5f0 *
                                   metric_factor[column, NNODES, layer]
    end

    if layer > 1
        theta_m = u_vertical[column, 3, NNODES, layer - 1] /
                  u_vertical[column, 1, NNODES, layer - 1]
        dpdrhotheta_m = dp_factor *
                        u_vertical[column, 3, NNODES, layer - 1]^dp_exponent
        # Interface potential temperature of the dissipative flux of `rho_theta`
        theta_interface = 0.5f0 * (theta_m + theta_element[1])
        # Factors of the last node of the layer below (rows `index_boundary - 2` and
        # `index_boundary - 1`) and of the first node of this layer (rows `index_boundary`
        # and `index_boundary + 1`)
        fM = 0.5f0 * metric_factor[column, NNODES, layer - 1]
        f1 = 0.5f0 * metric_factor[column, 1, layer]
        # (index_boundary-1,index_boundary-1)
        schur_band[column, 4, index_boundary - 1] += theta_interface * dpdrhotheta_m *
                                                     inv_speed_of_sound *
                                                     inv_boundary_weight * fM
        # (index_boundary,index_boundary-1)
        schur_band[column, 5, index_boundary - 1] += -theta_interface * dpdrhotheta_m *
                                                     inv_speed_of_sound * f1 *
                                                     inv_boundary_weight
        # (index_boundary-1,index_boundary)
        schur_band[column, 3, index_boundary] += -theta_interface *
                                                 dpdrhotheta_element[1] *
                                                 inv_speed_of_sound *
                                                 inv_boundary_weight * fM
        # (index_boundary,index_boundary)
        schur_band[column, 4, index_boundary] += theta_interface *
                                                 dpdrhotheta_element[1] *
                                                 inv_speed_of_sound * f1 *
                                                 inv_boundary_weight
        # (index_boundary-2,index_boundary-2)
        schur_band[column, 4, index_boundary - 2] += speed_of_sound * fM *
                                                     inv_boundary_weight
        # (index_boundary+1,index_boundary-2)
        schur_band[column, 7, index_boundary - 2] += -speed_of_sound * f1 *
                                                     inv_boundary_weight
        # (index_boundary-2,index_boundary+1)
        schur_band[column, 1, index_boundary + 1] += -speed_of_sound *
                                                     inv_boundary_weight * fM
        # (index_boundary+1,index_boundary+1)
        schur_band[column, 4, index_boundary + 1] += speed_of_sound * f1 *
                                                     inv_boundary_weight
        # (index_boundary,index_boundary-2)
        schur_band[column, 6, index_boundary - 2] += -theta_m * f1 * inv_boundary_weight
        # (index_boundary+1,index_boundary-1)
        schur_band[column, 6, index_boundary - 1] += -dpdrhotheta_m * f1 *
                                                     inv_boundary_weight
        # (index_boundary-1,index_boundary+1)
        schur_band[column, 2, index_boundary + 1] += theta_element[1] *
                                                     inv_boundary_weight * fM
        # (index_boundary-2,index_boundary)
        schur_band[column, 2, index_boundary] += dpdrhotheta_element[1] *
                                                 inv_boundary_weight * fM
    end

    if layer == 1
        schur_band[column, 3, index_boundary + 1] += metric_factor[column,
                                                                   1,
                                                                   layer] *
                                                     derivative_hat[1, 1] *
                                                     theta_element[1]
        schur_band[column, 5, index_boundary] += -metric_factor[column, 1,
                                                                layer] *
                                                 derivative_hat[1, 1] *
                                                 dpdrhotheta_element[1]
        schur_band[column, 4, index_boundary + 1] += metric_factor[column,
                                                                   1,
                                                                   layer] *
                                                     speed_of_sound *
                                                     inv_boundary_weight
    end
    schur_band[column, 2, index_boundary + 2] += metric_factor[column, 1,
                                                               layer] *
                                                 derivative_hat[1, NNODES] *
                                                 theta_element[NNODES]
    schur_band[column, 2, index_boundary + 3] += metric_factor[column, 1,
                                                               layer] *
                                                 derivative_hat[1, NNODES] *
                                                 dpdrhotheta_element[NNODES]
    schur_band[column, 6, index_boundary + 1] += metric_factor[column,
                                                               NNODES,
                                                               layer] *
                                                 derivative_hat[NNODES, 1] *
                                                 theta_element[1]
    schur_band[column, 6, index_boundary] += metric_factor[column, NNODES,
                                                           layer] *
                                             derivative_hat[NNODES, 1] *
                                             dpdrhotheta_element[1]
    if layer == nlayers
        schur_band[column, 5, index_boundary + 2] += metric_factor[column,
                                                                   NNODES,
                                                                   layer] *
                                                     derivative_hat[NNODES, NNODES] *
                                                     theta_element[NNODES]
        schur_band[column, 3, index_boundary + 3] += -metric_factor[column,
                                                                    NNODES,
                                                                    layer] *
                                                     derivative_hat[NNODES, NNODES] *
                                                     dpdrhotheta_element[NNODES]
        schur_band[column, 4, index_boundary + 2] += metric_factor[column,
                                                                   NNODES,
                                                                   layer] *
                                                     speed_of_sound *
                                                     inv_boundary_weight
    end
    return nothing
end

@kernel inbounds=true function fill_jacobian_KAkernel!(nlayers, B1m_34, B1_1,
                                                       B1_4, B1p_12,
                                                       schur_element, schur_band,
                                                       @Const(u_vertical),
                                                       @Const(metric_factor),
                                                       @Const(derivative_hat),
                                                       speed_of_sound,
                                                       inv_speed_of_sound,
                                                       inv_boundary_weight, dp_factor,
                                                       dp_exponent, inv_gamma_dt,
                                                       gamma_dt,
                                                       @Const(gravity), theta,
                                                       dpdrhotheta,
                                                       ::Val{NNODES}) where {NNODES}
    column = @index(Global)
    theta_element = MVector{NNODES, eltype(u_vertical)}(undef)
    dpdrhotheta_element = MVector{NNODES, eltype(u_vertical)}(undef)
    schur_correction = MMatrix{NNODES, NNODES, eltype(u_vertical)}(undef)

    for layer in 1:nlayers
        fill_jacobian_element!(layer, column, nlayers, Val(NNODES), theta_element,
                               dpdrhotheta_element, schur_correction,
                               B1m_34, B1_1, B1_4, B1p_12,
                               schur_element, schur_band, u_vertical,
                               metric_factor,
                               derivative_hat,
                               speed_of_sound, inv_speed_of_sound,
                               inv_boundary_weight,
                               dp_factor, dp_exponent, inv_gamma_dt, gamma_dt,
                               gravity,
                               theta, dpdrhotheta)
    end
end

function fill_jacobian!(backend::Backend, jacobian::VerticalJacobianTwoSidedCache,
                        equations,
                        gamma_dt)
    (; B1m_34, B1_1, B1_4, B1p_12, schur_element, schur_band, gravity, theta,
    dpdrhotheta, ncolumns, metric_factor, derivative_hat, inv_boundary_weight,
    u_vertical) = jacobian
    speed_of_sound = jacobian.speed_of_sound
    nnodes_ = number_of_nodes(jacobian)
    nlayers = number_of_layers(jacobian)
    dp_factor = equations.gamma * equations.K
    dp_exponent = equations.gamma - 1
    inv_speed_of_sound = inv(speed_of_sound)
    inv_gamma_dt = inv(gamma_dt)
    # The fill and the boundary Schur complement accumulate into `schur_band`
    fill!(schur_band, zero(eltype(schur_band)))
    @views schur_band[:, 4, :] .= inv_gamma_dt

    kernel! = fill_jacobian_KAkernel!(backend)
    kernel!(nlayers, B1m_34, B1_1, B1_4, B1p_12, schur_element, schur_band, u_vertical,
            metric_factor,
            derivative_hat, speed_of_sound, inv_speed_of_sound,
            inv_boundary_weight, dp_factor, dp_exponent, inv_gamma_dt, gamma_dt,
            gravity,
            theta, dpdrhotheta, Val(nnodes_);
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Boundary Schur complement: the contributions of the interior nodes of every element,
# eliminated with the local Schur complements `schur_element`, to the band matrix
# `schur_band`

# Contribution `-C A^{-1} B` of element `layer` of the column `column`; the layers of a
# column are processed in order by the same thread. The two columns of `B` of the first
# and the last node are processed together, first those of `rho_theta` and then those of
# `rho_w`. With the blocks `A_ij` of the element block `A_e` of the reference in the
# docstring of `VerticalJacobianTwoSided`, `x = A^{-1} b` of a column `b` of `B` follows
# from the local Schur complement `S_e`:
#   S_e x_rho_w = b_rho_w - A31 A11^{-1} b_rho - A32 A22^{-1} b_rho_theta,
#   x_rho = A11^{-1} (b_rho - A13 x_rho_w),
#   x_rho_theta = A22^{-1} (b_rho_theta - A23 x_rho_w).
# The arrays `x_rho`, `x_rho_theta`, and `x_rho_w` first hold the components of `b` and
# are then overwritten by those of `x`; `x_rho_1` and `x_rho_M` hold the density in the
# first and the last node, and `schur_update` the resulting entries `-C x`.
Base.@propagate_inbounds function schur_complement_boundary_element!(layer, column,
                                                                     nlayers,
                                                                     ::Val{NNODES},
                                                                     x_rho,
                                                                     x_rho_theta,
                                                                     x_rho_w, x_rho_1,
                                                                     x_rho_M,
                                                                     schur_update,
                                                                     B1m_34,
                                                                     B1p_12, B1_4, B1_1,
                                                                     schur_element,
                                                                     schur_band,
                                                                     gamma_dt,
                                                                     gravity, theta,
                                                                     dpdrhotheta,
                                                                     metric_factor,
                                                                     derivative_hat) where {NNODES}
    RealT = eltype(schur_band)
    zero_ = zero(RealT)
    # Offset of the 4 boundary unknowns of this layer in the condensed system
    offset_boundary = (layer - 1) * 4
    # Columns of `rho_theta` of the first and the last node: `b_rho_w = A32[:, 1]` and
    # `A32[:, M]`, `b_rho_theta = 0`, and the density entries `B1_1` and `B1_4` in the edge
    # nodes, where `A31` vanishes; thus, `S_e x_rho_w = b_rho_w`
    for i in 1:(NNODES - 2)
        x_rho_w[i, 1] = (metric_factor[column, i + 1, layer] *
                         derivative_hat[i + 1, 1] *
                         dpdrhotheta[column, 1, layer])
        x_rho_w[i, 2] = (metric_factor[column, i + 1, layer] *
                         derivative_hat[i + 1, NNODES] *
                         dpdrhotheta[column, NNODES, layer])
    end

    # Solve with the LU factors of `S_e` in place: forward substitution
    for k in 1:(NNODES - 3)
        for i in (k + 1):(NNODES - 2)
            x_rho_w[i, 1] -= schur_element[column, i, k, layer] * x_rho_w[k, 1]
            x_rho_w[i, 2] -= schur_element[column, i, k, layer] * x_rho_w[k, 2]
        end
    end
    # Backward substitution
    for k in (NNODES - 2):-1:1
        x_rho_w[k, 1] /= schur_element[column, k, k, layer]
        x_rho_w[k, 2] /= schur_element[column, k, k, layer]
        for i in 1:(k - 1)
            x_rho_w[i, 1] -= schur_element[column, i, k, layer] * x_rho_w[k, 1]
            x_rho_w[i, 2] -= schur_element[column, i, k, layer] * x_rho_w[k, 2]
        end
    end

    # Density in the edge nodes
    x_rho_1[1] = -B1_1[column, layer]
    x_rho_1[2] = zero_
    x_rho_M[1] = zero_
    x_rho_M[2] = -B1_4[column, layer]
    for j in 1:(NNODES - 2)
        x_rho_1[1] += (metric_factor[column, 1, layer] *
                       derivative_hat[1, j + 1]) * x_rho_w[j, 1]
        x_rho_1[2] += (metric_factor[column, 1, layer] *
                       derivative_hat[1, j + 1]) * x_rho_w[j, 2]
        x_rho_M[1] += (metric_factor[column, NNODES, layer] *
                       derivative_hat[NNODES, j + 1]) *
                      x_rho_w[j, 1]
        x_rho_M[2] += (metric_factor[column, NNODES, layer] *
                       derivative_hat[NNODES, j + 1]) *
                      x_rho_w[j, 2]
    end
    x_rho_1[1] *= -gamma_dt
    x_rho_1[2] *= -gamma_dt
    x_rho_M[1] *= -gamma_dt
    x_rho_M[2] *= -gamma_dt

    # Interior potential temperature density with `b_rho_theta = 0`
    for i in 1:(NNODES - 2)
        x_rho_theta[i, 1] = zero_
        x_rho_theta[i, 2] = zero_
        for j in 1:(NNODES - 2)
            x_rho_theta[i, 1] += (metric_factor[column, i + 1, layer] *
                                  derivative_hat[i + 1, j + 1] *
                                  theta[column, j + 1, layer]) * x_rho_w[j, 1]
            x_rho_theta[i, 2] += (metric_factor[column, i + 1, layer] *
                                  derivative_hat[i + 1, j + 1] *
                                  theta[column, j + 1, layer]) * x_rho_w[j, 2]
        end
        x_rho_theta[i, 1] *= -gamma_dt
        x_rho_theta[i, 2] *= -gamma_dt
    end

    # Entries `-C x` in the rows of the boundary unknowns
    schur_update[1, 1] = zero_
    schur_update[2, 1] = zero_
    schur_update[3, 1] = zero_
    schur_update[4, 1] = zero_
    schur_update[1, 2] = zero_
    schur_update[2, 2] = zero_
    schur_update[3, 2] = zero_
    schur_update[4, 2] = zero_
    for j in 1:(NNODES - 2)
        schur_update[1, 1] += -(metric_factor[column, 1, layer] *
                                derivative_hat[1, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j, 1]
        schur_update[1, 2] += -(metric_factor[column, 1, layer] *
                                derivative_hat[1, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j, 2]
        schur_update[4, 1] += -(metric_factor[column, NNODES, layer] *
                                derivative_hat[NNODES, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j, 1]
        schur_update[4, 2] += -(metric_factor[column, NNODES, layer] *
                                derivative_hat[NNODES, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j, 2]
        schur_update[2, 1] += -(metric_factor[column, 1, layer] *
                                derivative_hat[1, j + 1] *
                                dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 1]
        schur_update[2, 2] += -(metric_factor[column, 1, layer] *
                                derivative_hat[1, j + 1] *
                                dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 2]
        schur_update[3, 1] += -(metric_factor[column, NNODES, layer] *
                                derivative_hat[NNODES, j + 1] *
                                dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 1]
        schur_update[3, 2] += -(metric_factor[column, NNODES, layer] *
                                derivative_hat[NNODES, j + 1] *
                                dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 2]
    end
    gravity_1 = gravity[column, 1, layer]
    gravity_M = gravity[column, NNODES, layer]
    schur_update[2, 1] = schur_update[2, 1] - gravity_1 * x_rho_1[1]
    schur_update[2, 2] = schur_update[2, 2] - gravity_1 * x_rho_1[2]
    schur_update[3, 1] = schur_update[3, 1] - gravity_M * x_rho_M[1]
    schur_update[3, 2] = schur_update[3, 2] - gravity_M * x_rho_M[2]
    for i in 1:4
        schur_band[column, i + 3, offset_boundary + 1] += schur_update[i, 1]
        schur_band[column, i, offset_boundary + 4] += schur_update[i, 2]
    end

    # Columns of `rho_w` of the first and the last node: `b_rho = A13[:, 1]` and
    # `A13[:, M]`, whose diagonal entries vanish at the interfaces between elements since
    # `S / 2` has a zero diagonal, `b_rho_theta = A23[:, 1]` and `A23[:, M]`, and
    # `b_rho_w = 0`
    for i in 1:NNODES
        x_rho[i, 1] = ((layer > 1 && i == 1) ? zero_ :
                       metric_factor[column, i, layer] * derivative_hat[i, 1])
        x_rho[i, 2] = ((layer < nlayers && i == NNODES) ? zero_ :
                       metric_factor[column, i, layer] *
                       derivative_hat[i, NNODES])
    end
    for i in 1:(NNODES - 2)
        x_rho_theta[i, 1] = (metric_factor[column, i + 1, layer] *
                             derivative_hat[i + 1, 1] *
                             theta[column, 1, layer])
        x_rho_theta[i, 2] = (metric_factor[column, i + 1, layer] *
                             derivative_hat[i + 1, NNODES] *
                             theta[column, NNODES, layer])
    end
    # Right-hand side of `S_e x_rho_w`
    for i in 1:(NNODES - 2)
        x_rho_w[i, 1] = gravity[column, i + 1, layer] * x_rho[i + 1, 1]
        x_rho_w[i, 2] = gravity[column, i + 1, layer] * x_rho[i + 1, 2]
        for j in 1:(NNODES - 2)
            x_rho_w[i, 1] += (metric_factor[column, i + 1, layer] *
                              derivative_hat[i + 1, j + 1] *
                              dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 1]
            x_rho_w[i, 2] += (metric_factor[column, i + 1, layer] *
                              derivative_hat[i + 1, j + 1] *
                              dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 2]
        end
        x_rho_w[i, 1] *= -gamma_dt
        x_rho_w[i, 2] *= -gamma_dt
    end

    # Solve with the LU factors of `S_e` in place: forward substitution
    for k in 1:(NNODES - 3)
        for i in (k + 1):(NNODES - 2)
            x_rho_w[i, 1] -= schur_element[column, i, k, layer] * x_rho_w[k, 1]
            x_rho_w[i, 2] -= schur_element[column, i, k, layer] * x_rho_w[k, 2]
        end
    end
    # Backward substitution
    for k in (NNODES - 2):-1:1
        x_rho_w[k, 1] /= schur_element[column, k, k, layer]
        x_rho_w[k, 2] /= schur_element[column, k, k, layer]
        for i in 1:(k - 1)
            x_rho_w[i, 1] -= schur_element[column, i, k, layer] * x_rho_w[k, 1]
            x_rho_w[i, 2] -= schur_element[column, i, k, layer] * x_rho_w[k, 2]
        end
    end

    # Density in the edge nodes
    x_rho_1[1] = x_rho[1, 1]
    x_rho_1[2] = x_rho[1, 2]
    x_rho_M[1] = x_rho[NNODES, 1]
    x_rho_M[2] = x_rho[NNODES, 2]
    for j in 1:(NNODES - 2)
        x_rho_1[1] += -(metric_factor[column, 1, layer] *
                        derivative_hat[1, j + 1]) * x_rho_w[j, 1]
        x_rho_1[2] += -(metric_factor[column, 1, layer] *
                        derivative_hat[1, j + 1]) * x_rho_w[j, 2]
        x_rho_M[1] += -(metric_factor[column, NNODES, layer] *
                        derivative_hat[NNODES, j + 1]) *
                      x_rho_w[j, 1]
        x_rho_M[2] += -(metric_factor[column, NNODES, layer] *
                        derivative_hat[NNODES, j + 1]) *
                      x_rho_w[j, 2]
    end
    x_rho_1[1] *= gamma_dt
    x_rho_1[2] *= gamma_dt
    x_rho_M[1] *= gamma_dt
    x_rho_M[2] *= gamma_dt

    # Interior potential temperature density
    for i in 1:(NNODES - 2)
        for j in 1:(NNODES - 2)
            x_rho_theta[i, 1] += -(metric_factor[column, i + 1, layer] *
                                   derivative_hat[i + 1, j + 1] *
                                   theta[column, j + 1, layer]) * x_rho_w[j, 1]
            x_rho_theta[i, 2] += -(metric_factor[column, i + 1, layer] *
                                   derivative_hat[i + 1, j + 1] *
                                   theta[column, j + 1, layer]) * x_rho_w[j, 2]
        end
        x_rho_theta[i, 1] *= gamma_dt
        x_rho_theta[i, 2] *= gamma_dt
    end
    schur_update[1, 1] = zero_
    schur_update[2, 1] = zero_
    schur_update[3, 1] = zero_
    schur_update[4, 1] = zero_
    schur_update[1, 2] = zero_
    schur_update[2, 2] = zero_
    schur_update[3, 2] = zero_
    schur_update[4, 2] = zero_
    for j in 1:(NNODES - 2)
        schur_update[1, 1] += -(metric_factor[column, 1, layer] *
                                derivative_hat[1, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j, 1]
        schur_update[1, 2] += -(metric_factor[column, 1, layer] *
                                derivative_hat[1, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j, 2]
        schur_update[4, 1] += -(metric_factor[column, NNODES, layer] *
                                derivative_hat[NNODES, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j, 1]
        schur_update[4, 2] += -(metric_factor[column, NNODES, layer] *
                                derivative_hat[NNODES, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j, 2]
        schur_update[2, 1] += -(metric_factor[column, 1, layer] *
                                derivative_hat[1, j + 1] *
                                dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 1]
        schur_update[2, 2] += -(metric_factor[column, 1, layer] *
                                derivative_hat[1, j + 1] *
                                dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 2]
        schur_update[3, 1] += -(metric_factor[column, NNODES, layer] *
                                derivative_hat[NNODES, j + 1] *
                                dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 1]
        schur_update[3, 2] += -(metric_factor[column, NNODES, layer] *
                                derivative_hat[NNODES, j + 1] *
                                dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j, 2]
    end

    schur_update[2, 1] = schur_update[2, 1] - gravity_1 * x_rho_1[1]
    schur_update[3, 1] = schur_update[3, 1] - gravity_M * x_rho_M[1]
    schur_update[2, 2] = schur_update[2, 2] - gravity_1 * x_rho_1[2]
    schur_update[3, 2] = schur_update[3, 2] - gravity_M * x_rho_M[2]
    for i in 1:4
        schur_band[column, i + 2, offset_boundary + 2] += schur_update[i, 1]
        schur_band[column, i + 1, offset_boundary + 3] += schur_update[i, 2]
    end
    # Gravity couples the density in the edge nodes to the boundary unknowns of the
    # neighboring elements through `B1m_34` and `B1p_12`
    if layer > 1
        schur_band[column, 7, offset_boundary - 1] += -gravity_1 * gamma_dt *
                                                      B1m_34[column, 1, layer]
        schur_band[column, 6, offset_boundary] += -gravity_1 * gamma_dt *
                                                  B1m_34[column, 2, layer]
    end
    if layer < nlayers
        schur_band[column, 2, offset_boundary + 5] += -gravity_M * gamma_dt *
                                                      B1p_12[column, 1, layer]
        schur_band[column, 1, offset_boundary + 6] += -gravity_M * gamma_dt *
                                                      B1p_12[column, 2, layer]
    end
    return nothing
end

@kernel inbounds=true function schur_complement_boundary_KAkernel!(nlayers,
                                                                   @Const(B1m_34),
                                                                   @Const(B1p_12),
                                                                   @Const(B1_4),
                                                                   @Const(B1_1),
                                                                   @Const(schur_element),
                                                                   schur_band,
                                                                   gamma_dt,
                                                                   @Const(gravity),
                                                                   @Const(theta),
                                                                   @Const(dpdrhotheta),
                                                                   @Const(metric_factor),
                                                                   @Const(derivative_hat),
                                                                   ::Val{NNODES}) where {NNODES}
    column = @index(Global)
    RealT = eltype(schur_band)
    x_rho = MMatrix{NNODES, 2, RealT}(undef)
    x_rho_theta = MMatrix{NNODES - 2, 2, RealT}(undef)
    x_rho_w = MMatrix{NNODES - 2, 2, RealT}(undef)
    x_rho_1 = MVector{2, RealT}(undef)
    x_rho_M = MVector{2, RealT}(undef)
    schur_update = MMatrix{4, 2, RealT}(undef)

    for layer in 1:nlayers
        schur_complement_boundary_element!(layer, column, nlayers, Val(NNODES),
                                           x_rho,
                                           x_rho_theta,
                                           x_rho_w, x_rho_1, x_rho_M, schur_update,
                                           B1m_34, B1p_12,
                                           B1_4, B1_1, schur_element, schur_band,
                                           gamma_dt, gravity, theta,
                                           dpdrhotheta,
                                           metric_factor,
                                           derivative_hat)
    end
end

function schur_complement_boundary!(backend::Backend,
                                    jacobian::VerticalJacobianTwoSidedCache,
                                    gamma_dt)
    (; B1m_34, B1p_12, B1_4, B1_1, schur_element, ncolumns, schur_band, gravity, theta,
    dpdrhotheta, metric_factor, derivative_hat) = jacobian
    nnodes_ = number_of_nodes(jacobian)
    nlayers = number_of_layers(jacobian)

    kernel! = schur_complement_boundary_KAkernel!(backend)
    kernel!(nlayers, B1m_34, B1p_12, B1_4, B1_1, schur_element, schur_band, gamma_dt,
            gravity,
            theta, dpdrhotheta, metric_factor, derivative_hat,
            Val(nnodes_);
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Banded LU factorization of the boundary Schur complement and the corresponding solve.
# `schur_band[column, :, :]` stores the 3 lower and 3 upper diagonals in the usual banded
# layout, i.e., the entry `(i, j)` at `schur_band[column, i - j + 4, j]`.

Base.@propagate_inbounds function lu_band_column!(column, A, l, u, ::Val{N}) where {N}
    up1 = u + 1
    up2 = u + 2
    uplp1 = u + l + 1
    for k in 1:(N - 1)
        for i in up2:uplp1
            A[column, i, k] /= A[column, up1, k]
        end
        for i in up2:uplp1
            for j in (k + 1):min(k + u, N)
                A[column, i + k - j, j] -= A[column, i, k] * A[column, k + up1 - j, j]
            end
        end
    end
    return nothing
end

@kernel inbounds=true function lu_band_KAkernel!(schur_band, ::Val{N}) where {N}
    column = @index(Global)

    lu_band_column!(column, schur_band, 3, 3, Val(N))
end

function lu_band!(backend::Backend, jacobian::VerticalJacobianTwoSidedCache)
    (; ncolumns, schur_band) = jacobian
    kernel! = lu_band_KAkernel!(backend)
    kernel!(schur_band, Val(4 * number_of_layers(jacobian));
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

# Solve with the banded LU factors `A = schur_band` and the right-hand side `b = schur_rhs`
Base.@propagate_inbounds function ldiv_band_column!(column, A, b, l, u,
                                                    ::Val{N}) where {N}
    # Ly = b
    up1 = u + 1
    uplp1 = up1 + l

    for k in 1:(N - 1)
        for i in (k + 1):min(k + l, N)
            b[column, i] -= A[column, i - k + up1, k] * b[column, k]
        end
    end
    # Ux = y
    for k in N:-1:1
        b[column, k] /= A[column, up1, k]
        for i in max(k - u, 1):(k - 1)
            b[column, i] -= A[column, up1 + i - k, k] * b[column, k]
        end
    end
    return nothing
end

@kernel inbounds=true function ldiv_band_KAkernel!(@Const(schur_band), schur_rhs,
                                                   ::Val{N}) where {N}
    column = @index(Global)

    ldiv_band_column!(column, schur_band, schur_rhs, 3, 3, Val(N))
end

function ldiv_band!(backend::Backend, jacobian::VerticalJacobianTwoSidedCache)
    (; ncolumns, schur_band, schur_rhs) = jacobian
    kernel! = ldiv_band_KAkernel!(backend)
    kernel!(schur_band, schur_rhs, Val(4 * number_of_layers(jacobian));
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Forward substitution: the right-hand side of the condensed system
#   `schur_rhs = b_D - C * A^{-1} * b_A`,
# where `A` is the block of the interior and `D` the block of the boundary degrees of
# freedom

# Entries of `schur_rhs` of element `layer` of the column `column`, following the
# implementation in CGDycore.jl by Oswald Knoth
# (https://github.com/CliMA/CGDycore.jl). `x_rho`, `x_rho_theta`, and `x_rho_w` hold the
# components of `A^{-1} b_A` of the three variables, `x_rho_1` and `x_rho_M` those of the
# density in the first and the last node, and `schur_rhs_element` the entries of
# `b_D - C A^{-1} b_A`.
Base.@propagate_inbounds function ldiv_interior_forward_element!(layer, column, nlayers,
                                                                 ::Val{NNODES},
                                                                 x_rho, x_rho_theta,
                                                                 x_rho_w,
                                                                 schur_rhs_element,
                                                                 schur_element,
                                                                 schur_rhs, b, gamma_dt,
                                                                 gravity, theta,
                                                                 dpdrhotheta,
                                                                 metric_factor,
                                                                 derivative_hat) where {NNODES}
    RealT = eltype(b)
    # Offset of the 4 boundary unknowns of this layer in the condensed system
    offset_boundary = (layer - 1) * 4
    schur_rhs_element[1] = b[column, 3, 1, layer]
    schur_rhs_element[2] = b[column, 2, 1, layer]
    schur_rhs_element[3] = b[column, 2, NNODES, layer]
    schur_rhs_element[4] = b[column, 3, NNODES, layer]
    for i in 1:NNODES
        x_rho[i] = b[column, 1, i, layer]
    end
    for i in 1:(NNODES - 2)
        x_rho_theta[i] = b[column, 3, i + 1, layer]
        x_rho_w[i] = b[column, 2, i + 1, layer]
    end
    for i in 1:(NNODES - 2)
        x_rho_w[i] += -gamma_dt * gravity[column, i + 1, layer] * x_rho[i + 1]
        for j in 1:(NNODES - 2)
            x_rho_w[i] += -gamma_dt *
                          (metric_factor[column, i + 1, layer] *
                           derivative_hat[i + 1, j + 1] *
                           dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j]
        end
    end
    ldiv_dense!(layer, column, schur_element, x_rho_w, Val(NNODES - 2))
    x_rho_1 = x_rho[1]
    x_rho_M = x_rho[NNODES]
    for j in 1:(NNODES - 2)
        x_rho_1 += -(metric_factor[column, 1, layer] *
                     derivative_hat[1, j + 1]) * x_rho_w[j]
        x_rho_M += -(metric_factor[column, NNODES, layer] *
                     derivative_hat[NNODES, j + 1]) *
                   x_rho_w[j]
    end
    x_rho_1 *= gamma_dt
    x_rho_M *= gamma_dt
    for i in 1:(NNODES - 2)
        for j in 1:(NNODES - 2)
            x_rho_theta[i] += -(metric_factor[column, i + 1, layer] *
                                derivative_hat[i + 1, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j]
        end
        x_rho_theta[i] *= gamma_dt
    end
    for j in 1:(NNODES - 2)
        theta_j = theta[column, j + 1, layer]
        dpdrhotheta_j = dpdrhotheta[column, j + 1, layer]
        schur_rhs_element[1] += -(metric_factor[column, 1, layer] *
                                  derivative_hat[1, j + 1] * theta_j) *
                                x_rho_w[j]
        schur_rhs_element[4] += -(metric_factor[column, NNODES, layer] *
                                  derivative_hat[NNODES, j + 1] *
                                  theta_j) * x_rho_w[j]
        schur_rhs_element[2] += -(metric_factor[column, 1, layer] *
                                  derivative_hat[1, j + 1] * dpdrhotheta_j) *
                                x_rho_theta[j]
        schur_rhs_element[3] += -(metric_factor[column, NNODES, layer] *
                                  derivative_hat[NNODES, j + 1] *
                                  dpdrhotheta_j) * x_rho_theta[j]
    end
    schur_rhs_element[2] += -gravity[column, 1, layer] * x_rho_1
    schur_rhs_element[3] += -gravity[column, NNODES, layer] * x_rho_M
    schur_rhs[column, offset_boundary + 1] = schur_rhs_element[1]
    schur_rhs[column, offset_boundary + 2] = schur_rhs_element[2]
    schur_rhs[column, offset_boundary + 3] = schur_rhs_element[3]
    schur_rhs[column, offset_boundary + 4] = schur_rhs_element[4]
    return nothing
end

@kernel inbounds=true function ldiv_interior_forward_KAkernel!(nlayers,
                                                               @Const(schur_element),
                                                               @Const(b), schur_rhs,
                                                               gamma_dt,
                                                               @Const(gravity),
                                                               @Const(theta),
                                                               @Const(dpdrhotheta),
                                                               @Const(metric_factor),
                                                               @Const(derivative_hat),
                                                               ::Val{NNODES}) where {NNODES}
    column, layer = @index(Global, NTuple)
    x_rho = MVector{NNODES, eltype(b)}(undef)
    x_rho_theta = MVector{NNODES - 2, eltype(b)}(undef)
    x_rho_w = MVector{NNODES - 2, eltype(b)}(undef)
    schur_rhs_element = MVector{4, eltype(b)}(undef)

    ldiv_interior_forward_element!(layer, column, nlayers, Val(NNODES), x_rho,
                                   x_rho_theta, x_rho_w,
                                   schur_rhs_element, schur_element,
                                   schur_rhs, b, gamma_dt, gravity, theta,
                                   dpdrhotheta, metric_factor,
                                   derivative_hat)
end

function ldiv_interior_forward!(backend::Backend,
                                jacobian::VerticalJacobianTwoSidedCache, b,
                                gamma_dt)
    (; schur_rhs, schur_element, ncolumns, theta, dpdrhotheta, metric_factor, derivative_hat,
    gravity) = jacobian
    nnodes_ = number_of_nodes(jacobian)
    nlayers = number_of_layers(jacobian)

    kernel! = ldiv_interior_forward_KAkernel!(backend)
    kernel!(nlayers, schur_element, b, schur_rhs, gamma_dt, gravity, theta, dpdrhotheta,
            metric_factor, derivative_hat, Val(nnodes_);
            ndrange = (ncolumns, nlayers), layer_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Backward substitution: the solution of the boundary degrees of freedom `schur_rhs` and
# of the interior degrees of freedom

# Solution of element `layer` of the column `column`; `x_rho`, `x_rho_theta`, and `x_rho_w`
# hold the interior components of the three variables, and `rho_theta_1`, `rho_w_1`,
# `rho_w_M`, and `rho_theta_M` the boundary components from `schur_rhs`
Base.@propagate_inbounds function ldiv_interior_backward_element!(layer, column,
                                                                  nlayers,
                                                                  ::Val{NNODES},
                                                                  x_rho, x_rho_theta,
                                                                  x_rho_w,
                                                                  B1m_34, B1_1, B1_4,
                                                                  B1p_12,
                                                                  schur_element,
                                                                  schur_rhs, b,
                                                                  gamma_dt, gravity,
                                                                  theta, dpdrhotheta,
                                                                  metric_factor,
                                                                  derivative_hat) where {NNODES}
    RealT = eltype(b)
    zero_ = zero(RealT)
    # Offset of the 4 boundary unknowns of this layer in the condensed system
    offset_boundary = (layer - 1) * 4
    theta_1 = theta[column, 1, layer]
    theta_M = theta[column, NNODES, layer]
    dpdrhotheta_1 = dpdrhotheta[column, 1, layer]
    dpdrhotheta_M = dpdrhotheta[column, NNODES, layer]
    b[column, 3, 1, layer] = schur_rhs[column, offset_boundary + 1]
    b[column, 2, 1, layer] = schur_rhs[column, offset_boundary + 2]
    b[column, 2, NNODES, layer] = schur_rhs[column, offset_boundary + 3]
    b[column, 3, NNODES, layer] = schur_rhs[column, offset_boundary + 4]
    for i in 1:NNODES
        x_rho[i] = b[column, 1, i, layer]
    end
    for i in 1:(NNODES - 2)
        x_rho_theta[i] = b[column, 3, i + 1, layer]
        x_rho_w[i] = b[column, 2, i + 1, layer]
    end
    x_rho[1] -= B1_1[column, layer] * schur_rhs[column, offset_boundary + 1]
    x_rho[NNODES] -= B1_4[column, layer] * schur_rhs[column, offset_boundary + 4]
    rho_theta_1 = schur_rhs[column, offset_boundary + 1]
    rho_w_1 = schur_rhs[column, offset_boundary + 2]
    rho_w_M = schur_rhs[column, offset_boundary + 3]
    rho_theta_M = schur_rhs[column, offset_boundary + 4]
    for i in 1:NNODES
        b1 = metric_factor[column, i, layer] * derivative_hat[i, 1]
        b2 = metric_factor[column, i, layer] * derivative_hat[i, NNODES]
        if layer > 1 && i == 1
            b1 = zero_
        end
        if layer < nlayers && i == NNODES
            b2 = zero_
        end
        x_rho[i] -= b1 * rho_w_1 + b2 * rho_w_M
    end

    for i in 1:(NNODES - 2)
        x_rho_theta[i] -= (metric_factor[column, i + 1, layer] *
                           derivative_hat[i + 1, 1] * theta_1) *
                          rho_w_1 +
                          (metric_factor[column, i + 1, layer] *
                           derivative_hat[i + 1, NNODES] *
                           theta_M) * rho_w_M
        x_rho_w[i] -= (metric_factor[column, i + 1, layer] *
                       derivative_hat[i + 1, 1] * dpdrhotheta_1) *
                      rho_theta_1 +
                      (metric_factor[column, i + 1, layer] *
                       derivative_hat[i + 1, NNODES] *
                       dpdrhotheta_M) * rho_theta_M
    end
    if layer > 1
        x_rho[1] -= B1m_34[column, 1, layer] * schur_rhs[column, offset_boundary - 1] +
                    B1m_34[column, 2, layer] * schur_rhs[column, offset_boundary]
    end
    if layer < nlayers
        x_rho[NNODES] -= B1p_12[column, 1, layer] *
                         schur_rhs[column, offset_boundary + 5] +
                         B1p_12[column, 2, layer] *
                         schur_rhs[column, offset_boundary + 6]
    end
    for i in 1:(NNODES - 2)
        x_rho_w[i] += -gamma_dt * gravity[column, i + 1, layer] * x_rho[i + 1]
        for j in 1:(NNODES - 2)
            x_rho_w[i] += -gamma_dt *
                          (metric_factor[column, i + 1, layer] *
                           derivative_hat[i + 1, j + 1] *
                           dpdrhotheta[column, j + 1, layer]) * x_rho_theta[j]
        end
    end
    ldiv_dense!(layer, column, schur_element, x_rho_w, Val(NNODES - 2))
    for i in 1:NNODES
        for j in 1:(NNODES - 2)
            x_rho[i] += -(metric_factor[column, i, layer] *
                          derivative_hat[i, j + 1]) *
                        x_rho_w[j]
        end
        x_rho[i] *= gamma_dt
    end
    for i in 1:(NNODES - 2)
        for j in 1:(NNODES - 2)
            x_rho_theta[i] += -(metric_factor[column, i + 1, layer] *
                                derivative_hat[i + 1, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j]
        end
        x_rho_theta[i] *= gamma_dt
    end
    for i in 1:NNODES
        b[column, 1, i, layer] = x_rho[i]
    end
    for i in 1:(NNODES - 2)
        b[column, 3, i + 1, layer] = x_rho_theta[i]
        b[column, 2, i + 1, layer] = x_rho_w[i]
    end
    return nothing
end

@kernel inbounds=true function ldiv_interior_backward_KAkernel!(nlayers,
                                                                @Const(B1m_34),
                                                                @Const(B1_1),
                                                                @Const(B1_4),
                                                                @Const(B1p_12),
                                                                @Const(schur_element),
                                                                @Const(schur_rhs),
                                                                b,
                                                                gamma_dt,
                                                                @Const(gravity),
                                                                @Const(theta),
                                                                @Const(dpdrhotheta),
                                                                @Const(metric_factor),
                                                                @Const(derivative_hat),
                                                                ::Val{NNODES}) where {NNODES}
    column, layer = @index(Global, NTuple)
    x_rho = MVector{NNODES, eltype(b)}(undef)
    x_rho_theta = MVector{NNODES - 2, eltype(b)}(undef)
    x_rho_w = MVector{NNODES - 2, eltype(b)}(undef)

    ldiv_interior_backward_element!(layer, column, nlayers, Val(NNODES), x_rho,
                                    x_rho_theta, x_rho_w,
                                    B1m_34,
                                    B1_1, B1_4, B1p_12, schur_element, schur_rhs, b,
                                    gamma_dt,
                                    gravity, theta, dpdrhotheta,
                                    metric_factor,
                                    derivative_hat)
end

function ldiv_interior_backward!(backend::Backend,
                                 jacobian::VerticalJacobianTwoSidedCache,
                                 b,
                                 gamma_dt)
    (; B1m_34, B1_1, B1_4, B1p_12, schur_element, schur_rhs, ncolumns, gravity, theta,
    dpdrhotheta, metric_factor, derivative_hat) = jacobian
    nnodes_ = number_of_nodes(jacobian)
    nlayers = number_of_layers(jacobian)

    kernel! = ldiv_interior_backward_KAkernel!(backend)
    kernel!(nlayers, B1m_34, B1_1, B1_4, B1p_12, schur_element, schur_rhs, b, gamma_dt,
            gravity,
            theta,
            dpdrhotheta, metric_factor, derivative_hat, Val(nnodes_);
            ndrange = (ncolumns, nlayers), layer_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Dense LU factorization and solve of the local Schur complement `S_e` of one element
# (indexed `[column, i, j, layer]`, `N = nnodes - 2`)
Base.@propagate_inbounds function lu_dense!(layer, column, A, ::Val{N}) where {N}
    for j in 1:(N - 1)
        # One division and `N - j` multiplications instead of `N - j` divisions
        inv_jj = inv(A[column, j, j, layer])
        for i in (j + 1):N
            A[column, i, j, layer] *= inv_jj
        end
        for k in (j + 1):N
            a_jk = A[column, j, k, layer]
            @simd for i in (j + 1):N
                A[column, i, k, layer] -= A[column, i, j, layer] * a_jk
            end
        end
    end
    return nothing
end

Base.@propagate_inbounds function ldiv_dense!(layer, column, A, b, ::Val{N}) where {N}
    # Forward loop
    for k in 1:(N - 1)
        for i in (k + 1):N
            b[i] -= A[column, i, k, layer] * b[k]
        end
    end
    # Backward loop
    for k in N:-1:1
        b[k] /= A[column, k, k, layer]
        for i in 1:(k - 1)
            b[i] -= A[column, i, k, layer] * b[k]
        end
    end
    return nothing
end
end # @muladd
