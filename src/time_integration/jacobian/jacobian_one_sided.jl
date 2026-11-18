# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Vertically implicit Jacobian with the one-sided condensation of the nested
# Schur-complement solver. The nodes `1, ..., M-1` of each element are condensed by a local
# Schur complement onto the vertical momentum (dense LU factors `schur_element`), and the
# last node of each element forms a block-tridiagonal matrix with 3x3 blocks
# (`schur_diagonal`, `schur_lower`, `schur_upper`), factorized by a block LU.
# The vertical operator is the same as the one of `JacobianVerticalTwoSided`, see the
# comment at the top of `jacobian_two_sided.jl`.

"""
    JacobianVerticalOneSided(; c = 340)

Approximate Jacobian of the vertical dynamics (vertically propagating acoustic waves and
gravity) of the [`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`](@ref) on a
[`P4estMeshCubedSphereTopography`](@ref), used as linear operator of linearly implicit time
integration methods, i.e., Rosenbrock-W methods such as [`SSPRW32`](@ref) and linear IMEX
methods such as [`ARS343`](@ref), for horizontally explicit, vertically implicit (HEVI)
time integration, including the nested Schur-complement solver with the one-sided
condensation (block-tridiagonal reduced system) of
- Marco Artiano, Andrés Rueda-Ramírez, Oswald Knoth, Hendrik Ranocha, Gregor Gassner (2026)
  GPU-Accelerated Nested Schur-Complement Solvers for Atmospheric Flow Simulations using
  Horizontally-Explicit, Vertically-Implicit Schemes and Discontinuous Galerkin Methods
  [arXiv: 2610.15775](https://arxiv.org/abs/2610.15775)

The approximate Jacobian is the same as the one of [`JacobianVerticalTwoSided`](@ref), see
there for its usage and the approximations of the dissipation and the gravity. Only the
condensation of the linear systems differs, so both solve the same linear systems up to
round-off errors.

### Permutation
The unknowns of a column are the density `rho`, the vertical momentum `rho_w`, i.e., the
component of the momentum in the unit contravariant direction of the vertical reference
coordinate, and the potential temperature density `rho_theta`. For a column of `nz`
elements of polynomial degree `polydeg` (`M = polydeg + 1` LGL nodes), the unknowns are
grouped into a *local* block per element and a *global* (edge) block shared by the whole
column:

* the `rho` DOFs `1..M-1` of element 1, then
* the `rho_theta` DOFs `1..M-1` of element 1, then
* the `rho_w` DOFs `1..M-1` of element 1, then
* ... (repeated, element by element) ...
* finally, 3 edge DOFs per element, ordered element by element:
  `rho_M, rho_w_M, rho_theta_M`.

All three variables are condensed, but only at the *last* (upper) node of each element --
that is what "one-sided" refers to.

### Sparsity pattern
Example for `nz = 4`, `polydeg = 5` (`M = 6`), in the permuted ordering above:
```
⠑⢄⠀⠀⠀⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀
⠀⠀⠑⢄⠀⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀
⢄⠀⢠⣤⣵⢟⠛⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀
⠀⠑⠼⠿⠿⠀⠑⢄⠀⠀⡀⠀⢀⣀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡇⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⢸⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⢸⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠁⡇⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠐⢄⠀⣶⣶⡟⢍⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠂⢱⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠛⠛⠃⠀⠑⢄⠀⠠⠀⠀⣤⣤⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠼⢠⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⢸⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀ = ⎡A B⎤
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⢸⣿⣿⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⡇⠀   ⎣C D⎦
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠀⠀⠑⢄⠀⠂⠀⢰⣶⣶⠀⠀⠀⠐⠃⡆
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⢸⣿⣿⠀⠀⠀⠠⠄⡇
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⣀⣑⣼⠿⠿⠀⠀⠀⢀⡀⢇
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⣿⣿⡇⠑⢄⠀⠀⠀⠀⠀⢸
⠀⠀⠐⠒⠒⠭⠭⠅⠀⠀⠇⠀⢸⣀⣀⠀⠀⢀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠻⢇⣀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠑⠒⠒⠀⠀⢘⣀⣀⠧⠤⠄⠀⠀⡄⠀⢠⠀⠀⠀⠈⠛⣤⡄⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠁⠀⠀⠥⠤⢜⣒⣒⠀⠀⠀⠈⠱⣦
```

`A` (block-diagonal, one block per element) holds, for the nodes `1..M-1`, the diagonals
of `rho`, `rho_theta`, and `rho_w`, which are `I / gamma_dt` except for the first node at an
interface between elements, where the LMARS flux adds dissipation, the dense coupling (LGL
derivative operator) of `rho` and `rho_theta` to `rho_w` (mass and potential temperature
fluxes) and of `rho_w` to `rho_theta` (pressure gradient), the pointwise coupling of `rho_w`
to `rho` (gravity), and, at the first node at an interface, the coupling of `rho` to
`rho_theta` through the LMARS mass flux.

`B` connects the local DOFs to the edge `rho_w_M`/`rho_theta_M` of the same element (dense
LGL coupling) and, through the LMARS flux at the lower interface, the first node of the
element to the edge DOFs of the element below. `C` connects the edge DOFs to the local
`rho_w`/`rho_theta` of the same element (dense LGL coupling) and, through the LMARS flux at
the upper interface, to the first node of the element above.

`D` is block-diagonal with 3x3 blocks: the edge DOFs of an element couple only to each
other (gravity, the LGL coupling at the wall, and the dissipation of the LMARS flux at the
upper interface), since the neighboring edge DOFs are separated by the first node of the
element above.

### Algorithm
The matrix `I / gamma_dt - J` is assembled and factorized by `update_jacobian!`, i.e.,
whenever the Jacobian is updated (see the keyword argument `update_interval` of
`Trixi.solve`), and the linear systems of the stages are solved by `solve_jacobian!`:
1. **Local Schur complement** (per element, in `update_jacobian!`): eliminate the
   diagonal blocks of `rho` and `rho_theta` analytically, leaving a dense `(M-1)x(M-1)`
   complement `S_e` on `rho_w`, which is LU-factorized without pivoting (stored in
   `schur_element`).
2. **Global Schur complement** (per column, in `update_jacobian!`): eliminating the
   local DOFs of every element with `S_e` leaves a `3*nz x 3*nz` block-tridiagonal system
   `S = D - C A^{-1} B` with 3x3 blocks on the edge DOFs. The first column of the
   off-diagonal blocks vanishes, since `rho_M` couples only to the edge DOFs of its own
   element. Thus, `rho_M` is eliminated within each diagonal block, and the remaining
   block-tridiagonal system with 2x2 blocks is factorized with a block LU without pivoting.
3. **Solve** (per stage, in `solve_jacobian!`): copy the right-hand side to the columns,
   eliminate the local DOFs with `S_e` (right-hand side `b_D - C A^{-1} b_A` of the
   condensed system), solve the block-tridiagonal system with the block LU factors of `S`
   for the edge DOFs, back-substitute into the local DOFs, and copy the solution back to
   Trixi.jl's layout. The horizontal components of the momentum are only multiplied by
   `gamma_dt`, i.e., solved with `J = 0`.

All kernels use KernelAbstractions.jl with the backend of the arrays of the
semidiscretization, i.e., the CPU backend of KernelAbstractions.jl for `Array`s
(independently of the threading backend of Trixi.jl) or a GPU backend. Columns are
independent, so the sequential recurrences (global Schur complement, block LU
factorization and solve) use one thread per column, the copies between Trixi.jl's layout
and the columns use one thread per node, and the independent work per element (fill,
forward elimination, back-substitution) uses one thread per `(column, element)`. The
contributions of an element to the right-hand side of the condensed system of the element
below are added by a separate kernel, so that no atomic operations are needed.

!!! note
    The mesh must not be refined, since the elements of a column must be consecutive
    with the lowest layer first. For parallel meshes, the columns must be kept on one rank
    (`keep_columns_together = true` in [`P4estMeshCubedSphereTopography`](@ref)).
"""
struct JacobianVerticalOneSided{RealT <: Real}
    speed_of_sound::RealT
end

JacobianVerticalOneSided(; c = 340) = JacobianVerticalOneSided(c)

# Arrays of the `JacobianVerticalOneSided` for a semidiscretization, created by
# `create_cache_jacobian`
struct JacobianVerticalOneSidedCache{NLAYERS, NNODES, RealT <: Real, BackendType,
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
    # and norm `n = |Ja^3|` (slot 4), laid out `[slot, i, j, k, element]`
    vertical_direction_and_norm::ArrayRealT5
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
    dp_drho_theta::ArrayRealT3
    # LU factors of the local Schur complements `S_e` of the elements in the nodes
    # `1..M-1`, `[column, node, node, layer]`
    schur_element::ArrayRealT4
    # Global Schur complement `S`, a block-tridiagonal matrix with the diagonal blocks
    # `[column, row, column_block, layer]` and the blocks below and above the diagonal
    # without their first, vanishing, column, `[column, row, column_block - 1, layer]`, where
    # `schur_lower[:, :, :, layer]` is the block `(layer + 1, layer)` and
    # `schur_upper[:, :, :, layer]` the block `(layer, layer + 1)`
    schur_diagonal::ArrayRealT4
    schur_lower::ArrayRealT4
    schur_upper::ArrayRealT4
    # Block LU factors of `S` (see `lu_block_per_column!`): the multipliers and the inverse
    # pivots of the elimination of `rho_M`, `[column, layer]`, and the 2x2 blocks of the
    # inverse diagonal, of the lower factor, and of the upper factor,
    # `[column, row, column_block, layer]`
    lu_multiplier::MatrixRealT
    lu_inv_pivot::MatrixRealT
    lu_inv_diagonal::ArrayRealT4
    lu_lower::ArrayRealT4
    lu_upper::ArrayRealT4
    # Right-hand side and solution of the condensed system, `[column, row]`, and the
    # contributions of every element to the right-hand side of the edge DOFs of the element
    # below, `[column, v, layer]`
    schur_rhs::MatrixRealT
    schur_rhs_lower::ArrayRealT3
end

@inline function nlayers(::JacobianVerticalOneSidedCache{NLAYERS}) where {NLAYERS}
    return NLAYERS
end

@inline function Trixi.nnodes(::JacobianVerticalOneSidedCache{<:Any, NNODES}) where {NNODES}
    return NNODES
end

# Check whether the arrays of the columns in `jacobian` have the axes we assume they must
# have in the kernels of the vertical Jacobian.
function check_axes_columns(jacobian::JacobianVerticalOneSidedCache)
    (; ncolumns, metric_factor, gravity, derivative_hat, u_vertical, theta,
    dp_drho_theta, schur_element, schur_diagonal, schur_lower, schur_upper, lu_multiplier,
    lu_inv_pivot, lu_inv_diagonal, lu_lower, lu_upper, schur_rhs,
    schur_rhs_lower) = jacobian
    columns = Base.OneTo(ncolumns)
    nodes = Base.OneTo(nnodes(jacobian))
    local_nodes = Base.OneTo(nnodes(jacobian) - 1)
    layers = Base.OneTo(nlayers(jacobian))
    interfaces = Base.OneTo(nlayers(jacobian) - 1)
    Trixi.check_axes(metric_factor, (columns, nodes, layers))
    Trixi.check_axes(gravity, (columns, nodes, layers))
    Trixi.check_axes(theta, (columns, nodes, layers))
    Trixi.check_axes(dp_drho_theta, (columns, nodes, layers))
    Trixi.check_axes(derivative_hat, (nodes, nodes))
    Trixi.check_axes(u_vertical, (columns, Base.OneTo(3), nodes, layers))
    Trixi.check_axes(schur_element, (columns, local_nodes, local_nodes, layers))
    Trixi.check_axes(schur_diagonal, (columns, Base.OneTo(3), Base.OneTo(3), layers))
    Trixi.check_axes(schur_lower, (columns, Base.OneTo(3), Base.OneTo(2), interfaces))
    Trixi.check_axes(schur_upper, (columns, Base.OneTo(3), Base.OneTo(2), interfaces))
    Trixi.check_axes(lu_multiplier, (columns, layers))
    Trixi.check_axes(lu_inv_pivot, (columns, layers))
    Trixi.check_axes(lu_inv_diagonal, (columns, Base.OneTo(2), Base.OneTo(2), layers))
    Trixi.check_axes(lu_lower, (columns, Base.OneTo(2), Base.OneTo(2), interfaces))
    Trixi.check_axes(lu_upper, (columns, Base.OneTo(2), Base.OneTo(2), interfaces))
    # 3 edge unknowns per layer
    Trixi.check_axes(schur_rhs, (columns, Base.OneTo(3 * length(layers))))
    Trixi.check_axes(schur_rhs_lower, (columns, Base.OneTo(3), layers))
    return nothing
end

function create_cache_jacobian(jacobian::JacobianVerticalOneSided, semi, u_ode)
    (; mesh, equations, solver, cache) = semi
    if !(equations isa CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
        throw(ArgumentError("`JacobianVerticalOneSided` requires the " *
                            "`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`"))
    end
    nlayers_ = nlayers(mesh)
    nnodes_ = nnodes(solver)
    nelements_ = nelements(solver, cache)
    # Compare the global numbers of elements and trees, since `nelements` is local to the
    # MPI rank for parallel meshes
    if Trixi.ncellsglobal(mesh) != Trixi.ntrees(mesh)
        throw(ArgumentError("`JacobianVerticalOneSided` requires a mesh " *
                            "without refinement"))
    end
    ncolumns = nnodes_^2 * nelements_ ÷ nlayers_
    RealT = eltype(cache.elements.node_coordinates)

    # The metric terms are computed on the CPU and moved to the storage type of the
    # semidiscretization
    storage = Trixi.storage_type(cache.elements)
    check_vertical_ordering(nlayers_, semi)
    index_maps = compute_element_index_maps(nlayers_, solver, cache)
    horizontal_index_of_element = Adapt.adapt(storage, index_maps[1])
    vertical_index_of_element = Adapt.adapt(storage, index_maps[2])
    vertical_direction_and_norm = Adapt.adapt(storage,
                                              compute_vertical_direction_and_norm(semi))
    metric_factor = compute_metric_factor(ncolumns, nlayers_, semi, index_maps...)
    metric_factor = Adapt.adapt(storage, metric_factor)
    gravity = compute_gravity(ncolumns, nlayers_, semi, index_maps...)
    gravity = Adapt.adapt(storage, gravity)
    inv_boundary_weight = solver.basis.inverse_weights[1]
    derivative_hat = solver.basis.derivative_hat

    # The remaining arrays are allocated with the backend of the semidiscretization
    backend = get_backend(cache.elements.node_coordinates)
    u_vertical = KernelAbstractions.zeros(backend, RealT, ncolumns, 3, nnodes_,
                                          nlayers_)
    theta = KernelAbstractions.zeros(backend, RealT, ncolumns, nnodes_, nlayers_)
    dp_drho_theta = KernelAbstractions.zeros(backend, RealT, ncolumns, nnodes_,
                                             nlayers_)
    schur_element = KernelAbstractions.zeros(backend, RealT, ncolumns, nnodes_ - 1,
                                             nnodes_ - 1, nlayers_)
    schur_diagonal = KernelAbstractions.zeros(backend, RealT, ncolumns, 3, 3, nlayers_)
    schur_lower = KernelAbstractions.zeros(backend, RealT, ncolumns, 3, 2, nlayers_ - 1)
    schur_upper = KernelAbstractions.zeros(backend, RealT, ncolumns, 3, 2, nlayers_ - 1)
    lu_multiplier = KernelAbstractions.zeros(backend, RealT, ncolumns, nlayers_)
    lu_inv_pivot = KernelAbstractions.zeros(backend, RealT, ncolumns, nlayers_)
    lu_inv_diagonal = KernelAbstractions.zeros(backend, RealT, ncolumns, 2, 2, nlayers_)
    lu_lower = KernelAbstractions.zeros(backend, RealT, ncolumns, 2, 2, nlayers_ - 1)
    lu_upper = KernelAbstractions.zeros(backend, RealT, ncolumns, 2, 2, nlayers_ - 1)
    schur_rhs = KernelAbstractions.zeros(backend, RealT, ncolumns, 3 * nlayers_)
    schur_rhs_lower = KernelAbstractions.zeros(backend, RealT, ncolumns, 3, nlayers_)

    return JacobianVerticalOneSidedCache{nlayers_, nnodes_, RealT, typeof(backend),
                                         typeof(schur_rhs), typeof(gravity),
                                         typeof(schur_element),
                                         typeof(vertical_direction_and_norm),
                                         typeof(horizontal_index_of_element),
                                         typeof(derivative_hat)}(ncolumns,
                                                                 convert(RealT,
                                                                         jacobian.speed_of_sound),
                                                                 inv_boundary_weight,
                                                                 backend,
                                                                 vertical_direction_and_norm,
                                                                 metric_factor,
                                                                 gravity,
                                                                 horizontal_index_of_element,
                                                                 vertical_index_of_element,
                                                                 derivative_hat,
                                                                 u_vertical,
                                                                 theta, dp_drho_theta,
                                                                 schur_element,
                                                                 schur_diagonal,
                                                                 schur_lower,
                                                                 schur_upper,
                                                                 lu_multiplier,
                                                                 lu_inv_pivot,
                                                                 lu_inv_diagonal,
                                                                 lu_lower, lu_upper,
                                                                 schur_rhs,
                                                                 schur_rhs_lower)
end

# Assemble and factorize `I / gamma_dt - J` from the state `u_ode`, see `jacobian.jl`
function update_jacobian!(jacobian::JacobianVerticalOneSidedCache, u_ode, semi,
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
    @trixi_timeit_ext backend timer() "block LU" lu_block!(backend, jacobian)

    return nothing
end

# Solve `(I / gamma_dt - J) x = b_ode` in place, see `jacobian.jl`
function solve_jacobian!(b_ode, jacobian::JacobianVerticalOneSidedCache, semi, gamma_dt)
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

    @trixi_timeit_ext backend timer() "block solve" ldiv_block!(backend, jacobian)

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
# Entries of the derivative blocks of `A`, `B`, and `C`, which are computed from
# `derivative_hat`, `metric_factor`, `theta`, and `dp_drho_theta` where they are needed.
# The entries `(1, 1)` of `A13`, `A23`, and `A32` differ from `metric_factor * D`: at an
# interface between elements, the central part of the LMARS flux cancels the boundary term
# of the derivative, and at the wall the flux of `rho` and `rho_theta` vanishes, while the
# pressure doubles the boundary term of the momentum. Thus, the entries `(1, 1)` are
# multiplied by `factor_11_a13 = factor_11_a23 = 1` and `factor_11_a32 = -1` in the lowest
# layer and by zero otherwise.

# Entry `(i, j)` of `A13` (`factor_11 = factor_11_a13`), of `A23` (multiplied by
# `theta_j`), or of `A32` (`factor_11 = factor_11_a32`, multiplied by `dp_drho_theta_j`)
@inline function derivative_entry(derivative_hat, metric_factor_i, factor_11, i, j)
    entry = metric_factor_i * derivative_hat[i, j]
    return ifelse((i == 1) & (j == 1), factor_11 * entry, entry)
end

@inline function factors_11(layer)
    factor_11_a13 = layer == 1
    factor_11_a32 = ifelse(layer == 1, -1, 0)
    return factor_11_a13, factor_11_a32
end

# Inverse of the diagonal of `A22` in the node `i`
@inline function inv_a22(inv_a22_1, gamma_dt, i)
    return ifelse(i == 1, inv_a22_1, gamma_dt)
end

# Entries of the LMARS flux at the lower interface of a layer, with the last node of the
# layer below as `ll` and the first node of the layer as `rr`, where `factor_ll` and
# `factor_rr` are the factors `metric_factor / (2 w_1)` of the rows of these nodes:
# `A12[1, 1]` (`rho` and `rho_theta`), the dissipation of `A22[1, 1]` and `A33[1, 1]`, the
# entries of `B` in the rows `rho`, `rho_w`, and `rho_theta` of the first node and the
# columns `rho_w_M` and `rho_theta_M` of the layer below, and the entries of `C` in the rows
# `rho_M`, `rho_w_M`, and `rho_theta_M` of the layer below and the columns `rho_w` and
# `rho_theta` of the first node. The dissipative flux of `rho_theta` uses the interface
# potential temperature, i.e., the mean of both sides.
@inline function lower_interface_entries(theta_ll, dp_drho_theta_ll, theta_rr,
                                         dp_drho_theta_rr, factor_ll, factor_rr,
                                         speed_of_sound, inv_speed_of_sound)
    theta_interface = 0.5f0 * (theta_ll + theta_rr)
    a12_1 = dp_drho_theta_rr * inv_speed_of_sound * factor_rr
    a22_1_dissipation = theta_interface * dp_drho_theta_rr * inv_speed_of_sound *
                        factor_rr
    a33_1_dissipation = speed_of_sound * factor_rr
    b = SMatrix{3, 2}(-factor_rr,
                      -speed_of_sound * factor_rr,
                      -theta_ll * factor_rr,
                      -dp_drho_theta_ll * inv_speed_of_sound * factor_rr,
                      -dp_drho_theta_ll * factor_rr,
                      -theta_interface * dp_drho_theta_ll * inv_speed_of_sound *
                      factor_rr)
    c = SMatrix{3, 2}(factor_ll,
                      -speed_of_sound * factor_ll,
                      theta_rr * factor_ll,
                      -dp_drho_theta_rr * inv_speed_of_sound * factor_ll,
                      dp_drho_theta_rr * factor_ll,
                      -theta_interface * dp_drho_theta_rr * inv_speed_of_sound *
                      factor_ll)
    return a12_1, a22_1_dissipation, a33_1_dissipation, b, c
end

# `A12[1, 1]`, the inverse of `A22[1, 1]`, and the entries of `B` and `C` at the lower
# interface of element `layer` (see `lower_interface_entries`) from the potential
# temperature and `dp/d(rho_theta)` stored by the fill. At the wall, the entries of `B` and
# `C` vanish and `A12[1, 1]` and `A22[1, 1]` are those of the other nodes.
Base.@propagate_inbounds function first_node_entries(theta, dp_drho_theta,
                                                     metric_factor, speed_of_sound,
                                                     inv_speed_of_sound,
                                                     inv_boundary_weight, gamma_dt,
                                                     ::Val{NNODES}, column,
                                                     layer) where {NNODES}
    RealT = eltype(theta)
    if layer == 1
        zero_ = zero(SMatrix{3, 2, RealT})
        return zero(RealT), gamma_dt, zero_, zero_
    end
    factor_ll = 0.5f0 * metric_factor[column, NNODES, layer - 1] * inv_boundary_weight
    factor_rr = 0.5f0 * metric_factor[column, 1, layer] * inv_boundary_weight
    a12_1, a22_1, _, b, c = lower_interface_entries(theta[column, NNODES, layer - 1],
                                                    dp_drho_theta[column, NNODES,
                                                                  layer - 1],
                                                    theta[column, 1, layer],
                                                    dp_drho_theta[column, 1, layer],
                                                    factor_ll, factor_rr,
                                                    speed_of_sound, inv_speed_of_sound)
    return a12_1, inv(inv(gamma_dt) + a22_1), b, c
end

###############################################################################
# Fill of the Jacobian: the potential temperature and `dp/d(rho_theta)` per node, the LU
# factors `schur_element` of the local Schur complements, and the blocks of `D` in
# `schur_diagonal`. The blocks below and above the diagonal of `D` vanish. The entries of
# `A` in the first node and of `B` and `C` at the interfaces are computed where they are
# needed.

# Fill of element `layer` of the column `column`, which writes only the entries of this
# element and reads the state of the neighboring elements
Base.@propagate_inbounds function fill_jacobian_one_sided_per_element!(schur_element,
                                                                       schur_diagonal,
                                                                       schur_lower,
                                                                       schur_upper,
                                                                       theta,
                                                                       dp_drho_theta,
                                                                       theta_element,
                                                                       dp_drho_theta_element,
                                                                       u_vertical,
                                                                       metric_factor,
                                                                       gravity,
                                                                       derivative_hat,
                                                                       speed_of_sound,
                                                                       inv_speed_of_sound,
                                                                       inv_boundary_weight,
                                                                       dp_factor,
                                                                       dp_exponent,
                                                                       inv_gamma_dt,
                                                                       gamma_dt,
                                                                       ::Val{NLAYERS},
                                                                       ::Val{NNODES},
                                                                       column,
                                                                       layer) where {
                                                                                     NLAYERS,
                                                                                     NNODES
                                                                                     }
    RealT = eltype(u_vertical)
    zero_ = zero(RealT)
    for i in 1:NNODES
        theta_element[i] = u_vertical[column, 3, i, layer] /
                           u_vertical[column, 1, i, layer]
        dp_drho_theta_element[i] = dp_factor *
                                   u_vertical[column, 3, i, layer]^dp_exponent
        theta[column, i, layer] = theta_element[i]
        dp_drho_theta[column, i, layer] = dp_drho_theta_element[i]
    end

    # Entries of `A` in the first node, with the dissipation of the LMARS flux at the lower
    # interface or of the pressure at the wall. The state of the layer below is computed
    # from `u_vertical`, since another thread fills that layer.
    if layer == 1
        a12_1 = zero_
        inv_a22_1 = gamma_dt
        a33_1 = inv_gamma_dt +
                metric_factor[column, 1, layer] * speed_of_sound * inv_boundary_weight
    else
        theta_ll = u_vertical[column, 3, NNODES, layer - 1] /
                   u_vertical[column, 1, NNODES, layer - 1]
        dp_drho_theta_ll = dp_factor *
                           u_vertical[column, 3, NNODES, layer - 1]^dp_exponent
        factor_ll = 0.5f0 * metric_factor[column, NNODES, layer - 1] *
                    inv_boundary_weight
        factor_rr = 0.5f0 * metric_factor[column, 1, layer] * inv_boundary_weight
        a12_1, a22_1, a33_1, _, _ = lower_interface_entries(theta_ll, dp_drho_theta_ll,
                                                            theta_element[1],
                                                            dp_drho_theta_element[1],
                                                            factor_ll, factor_rr,
                                                            speed_of_sound,
                                                            inv_speed_of_sound)
        inv_a22_1 = inv(inv_gamma_dt + a22_1)
        a33_1 += inv_gamma_dt
    end

    # Local Schur complement
    #   S_e = A33 - A31 A11^{-1} A13 - (A32 - A31 A11^{-1} A12) A22^{-1} A23
    # with `A11 = I / gamma_dt`, `A31 = diag(gravity)`, and `A12` and the diagonal matrices
    # `A22` and `A33` as above
    factor_11_a13, factor_11_a32 = factors_11(layer)
    for j in 1:(NNODES - 1)
        for i in 1:(NNODES - 1)
            metric_factor_i = metric_factor[column, i, layer]
            # Entry `(A32 A22^{-1} A23)[i, j]`
            correction = zero_
            for k in 1:(NNODES - 1)
                a32 = derivative_entry(derivative_hat, metric_factor_i, factor_11_a32,
                                       i, k) * dp_drho_theta_element[k]
                a23 = derivative_entry(derivative_hat, metric_factor[column, k, layer],
                                       factor_11_a13, k, j) * theta_element[j]
                correction += a32 * a23 * inv_a22(inv_a22_1, gamma_dt, k)
            end
            a13 = derivative_entry(derivative_hat, metric_factor_i, factor_11_a13, i, j)
            a33 = ifelse(i == j, ifelse(i == 1, a33_1, inv_gamma_dt), zero_)
            schur_element[column, i, j, layer] = a33 -
                                                 gravity[column, i, layer] * gamma_dt *
                                                 a13 - correction
        end
        # `A31 A11^{-1} A12 A22^{-1} A23` has only the row of the first node
        a23_1j = derivative_entry(derivative_hat, metric_factor[column, 1, layer],
                                  factor_11_a13, 1, j) * theta_element[j]
        schur_element[column, 1, j, layer] += gravity[column, 1, layer] * gamma_dt *
                                              a12_1 * inv_a22_1 * a23_1j
    end
    lu_dense!(schur_element, Val(NNODES - 1), column, layer)

    # Block of `D` in the rows and columns `rho_M`, `rho_w_M`, and `rho_theta_M`
    metric_factor_M = metric_factor[column, NNODES, layer]
    theta_M = theta_element[NNODES]
    dp_drho_theta_M = dp_drho_theta_element[NNODES]
    d11 = inv_gamma_dt
    d21 = gravity[column, NNODES, layer]
    d31 = zero_
    if layer < NLAYERS
        # Dissipation of the LMARS flux at the upper interface, with this layer as `ll` and
        # the layer above as `rr`; the central part cancels the boundary term of the
        # derivative
        theta_rr = u_vertical[column, 3, 1, layer + 1] /
                   u_vertical[column, 1, 1, layer + 1]
        theta_interface = 0.5f0 * (theta_M + theta_rr)
        factor_ll = 0.5f0 * metric_factor_M * inv_boundary_weight
        d12 = zero_
        d22 = inv_gamma_dt + speed_of_sound * factor_ll
        d32 = zero_
        d13 = dp_drho_theta_M * inv_speed_of_sound * factor_ll
        d23 = zero_
        d33 = inv_gamma_dt +
              theta_interface * dp_drho_theta_M * inv_speed_of_sound *
              factor_ll
    else
        # Wall
        derivative_MM = metric_factor_M * derivative_hat[NNODES, NNODES]
        d12 = derivative_MM
        d22 = inv_gamma_dt + metric_factor_M * speed_of_sound * inv_boundary_weight
        d32 = derivative_MM * theta_M
        d13 = zero_
        d23 = -derivative_MM * dp_drho_theta_M
        d33 = inv_gamma_dt
    end
    schur_diagonal[column, 1, 1, layer] = d11
    schur_diagonal[column, 2, 1, layer] = d21
    schur_diagonal[column, 3, 1, layer] = d31
    schur_diagonal[column, 1, 2, layer] = d12
    schur_diagonal[column, 2, 2, layer] = d22
    schur_diagonal[column, 3, 2, layer] = d32
    schur_diagonal[column, 1, 3, layer] = d13
    schur_diagonal[column, 2, 3, layer] = d23
    schur_diagonal[column, 3, 3, layer] = d33

    # The boundary Schur complement accumulates into the blocks off the diagonal
    if layer < NLAYERS
        for k in 1:2, i in 1:3
            schur_lower[column, i, k, layer] = zero_
            schur_upper[column, i, k, layer] = zero_
        end
    end
    return nothing
end

@kernel inbounds=true function fill_jacobian_one_sided_KAkernel!(schur_element,
                                                                 schur_diagonal,
                                                                 schur_lower,
                                                                 schur_upper, theta,
                                                                 dp_drho_theta,
                                                                 @Const(u_vertical),
                                                                 @Const(metric_factor),
                                                                 @Const(gravity),
                                                                 @Const(derivative_hat),
                                                                 speed_of_sound,
                                                                 inv_speed_of_sound,
                                                                 inv_boundary_weight,
                                                                 dp_factor, dp_exponent,
                                                                 inv_gamma_dt, gamma_dt,
                                                                 ::Val{NLAYERS},
                                                                 ::Val{NNODES}) where {
                                                                                       NLAYERS,
                                                                                       NNODES
                                                                                       }
    column, layer = @index(Global, NTuple)
    theta_element = MVector{NNODES, eltype(u_vertical)}(undef)
    dp_drho_theta_element = MVector{NNODES, eltype(u_vertical)}(undef)

    fill_jacobian_one_sided_per_element!(schur_element, schur_diagonal, schur_lower,
                                         schur_upper, theta, dp_drho_theta,
                                         theta_element, dp_drho_theta_element,
                                         u_vertical,
                                         metric_factor, gravity, derivative_hat,
                                         speed_of_sound, inv_speed_of_sound,
                                         inv_boundary_weight, dp_factor, dp_exponent,
                                         inv_gamma_dt, gamma_dt, Val(NLAYERS),
                                         Val(NNODES), column, layer)
end

function fill_jacobian!(backend::Backend, jacobian::JacobianVerticalOneSidedCache,
                        equations, gamma_dt)
    (; schur_element, schur_diagonal, schur_lower, schur_upper, theta, dp_drho_theta,
    u_vertical, metric_factor, gravity, derivative_hat, inv_boundary_weight,
    ncolumns) = jacobian
    speed_of_sound = jacobian.speed_of_sound
    nnodes_ = nnodes(jacobian)
    nlayers_ = nlayers(jacobian)
    dp_factor = equations.gamma * equations.K
    dp_exponent = equations.gamma - 1
    inv_speed_of_sound = inv(speed_of_sound)
    inv_gamma_dt = inv(gamma_dt)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck check_axes_columns(jacobian)

    kernel! = fill_jacobian_one_sided_KAkernel!(backend)
    kernel!(schur_element, schur_diagonal, schur_lower, schur_upper, theta,
            dp_drho_theta, u_vertical, metric_factor, gravity, derivative_hat,
            speed_of_sound, inv_speed_of_sound, inv_boundary_weight, dp_factor,
            dp_exponent, inv_gamma_dt, gamma_dt, Val(nlayers_), Val(nnodes_);
            ndrange = (ncolumns, nlayers_), layer_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Boundary Schur complement: the contributions of the nodes `1..M-1` of every element,
# eliminated with the local Schur complements `schur_element`, to the block-tridiagonal
# matrix `schur_diagonal`, `schur_lower`, `schur_upper`

# Solve `A_e x = b` with the local Schur complement `S_e` for the right-hand sides `b` in the
# columns `r` of `x_rho_theta[:, r]` and `x_rho_w[:, r]`, which are solved together so that
# the entries of `A_e` and `S_e` are loaded once. With the blocks `A_ij` of `A_e`, `x`
# follows from
#   S_e x_rho_w = b_rho_w - A31 A11^{-1} b_rho - (A32 - A31 A11^{-1} A12) A22^{-1} b_rho_theta,
#   x_rho_theta = A22^{-1} (b_rho_theta - A23 x_rho_w).
# On entry, `x_rho_theta` holds `b_rho_theta` and `x_rho_w` holds
# `b_rho_w - A31 A11^{-1} b_rho`. Both are overwritten by the components of `x`. The density
# `x_rho` is not needed, since `C` does not couple to the density in the nodes `1..M-1`.
Base.@propagate_inbounds function eliminate_local_one_sided!(x_rho_theta, x_rho_w,
                                                             schur_element, gravity,
                                                             theta, dp_drho_theta,
                                                             metric_factor,
                                                             derivative_hat,
                                                             a12_1, inv_a22_1,
                                                             gamma_dt, ::Val{NNODES},
                                                             column,
                                                             layer) where {NNODES}
    factor_11_a13, factor_11_a32 = factors_11(layer)
    rhs = axes(x_rho_w, 2)
    # Right-hand side of `S_e x_rho_w`
    for j in 1:(NNODES - 1)
        inv_a22_j = inv_a22(inv_a22_1, gamma_dt, j)
        dp_drho_theta_j = dp_drho_theta[column, j, layer]
        for i in 1:(NNODES - 1)
            a32_inv_a22 = derivative_entry(derivative_hat,
                                           metric_factor[column, i, layer],
                                           factor_11_a32, i, j) *
                          dp_drho_theta_j * inv_a22_j
            for r in rhs
                x_rho_w[i, r] -= a32_inv_a22 * x_rho_theta[j, r]
            end
        end
    end
    a31_a12_inv_a22 = gravity[column, 1, layer] * gamma_dt * a12_1 * inv_a22_1
    for r in rhs
        x_rho_w[1, r] += a31_a12_inv_a22 * x_rho_theta[1, r]
    end

    # Solve `S_e x_rho_w` with the LU factors of `S_e` in place
    ldiv_dense_columns!(x_rho_w, schur_element, Val(NNODES - 1), column, layer)

    for j in 1:(NNODES - 1)
        theta_j = theta[column, j, layer]
        for i in 1:(NNODES - 1)
            a23 = derivative_entry(derivative_hat, metric_factor[column, i, layer],
                                   factor_11_a13, i, j) * theta_j
            for r in rhs
                x_rho_theta[i, r] -= a23 * x_rho_w[j, r]
            end
        end
    end
    for i in 1:(NNODES - 1)
        inv_a22_i = inv_a22(inv_a22_1, gamma_dt, i)
        for r in rhs
            x_rho_theta[i, r] *= inv_a22_i
        end
    end
    return nothing
end

# `ldiv_dense!` for the right-hand sides in the columns of `b`
Base.@propagate_inbounds function ldiv_dense_columns!(b, A, ::Val{N}, column,
                                                      layer) where {N}
    rhs = axes(b, 2)
    # Forward loop
    for k in 1:(N - 1)
        for i in (k + 1):N
            a_ik = A[column, i, k, layer]
            for r in rhs
                b[i, r] -= a_ik * b[k, r]
            end
        end
    end
    # Backward loop
    for k in N:-1:1
        a_kk = A[column, k, k, layer]
        for r in rhs
            b[k, r] /= a_kk
        end
        for i in 1:(k - 1)
            a_ik = A[column, i, k, layer]
            for r in rhs
                b[i, r] -= a_ik * b[k, r]
            end
        end
    end
    return nothing
end

# Entries `-C x` in the rows `rho_M`, `rho_w_M`, and `rho_theta_M` of element `layer`, with
# the column `r` of `x_rho_theta` and `x_rho_w` from `eliminate_local_one_sided!`
Base.@propagate_inbounds function minus_c_x_one_sided(x_rho_theta, x_rho_w, theta,
                                                      dp_drho_theta, metric_factor,
                                                      derivative_hat, ::Val{NNODES},
                                                      column, layer, r) where {NNODES}
    metric_factor_M = metric_factor[column, NNODES, layer]
    rho_M = rho_w_M = rho_theta_M = zero(eltype(x_rho_w))
    for j in 1:(NNODES - 1)
        derivative_Mj = metric_factor_M * derivative_hat[NNODES, j]
        rho_M -= derivative_Mj * x_rho_w[j, r]
        rho_w_M -= derivative_Mj * dp_drho_theta[column, j, layer] * x_rho_theta[j, r]
        rho_theta_M -= derivative_Mj * theta[column, j, layer] * x_rho_w[j, r]
    end
    return rho_M, rho_w_M, rho_theta_M
end

# Contribution `-C A^{-1} B` of element `layer` of the column `column`; the layers of a
# column are processed in order by the same thread, which adds the contributions of the
# columns of `B` of the edge DOFs of the layer below to the blocks of that layer. The
# column of `rho_M` of `B` vanishes.
Base.@propagate_inbounds function schur_complement_boundary_one_sided_per_element!(schur_diagonal,
                                                                                   schur_lower,
                                                                                   schur_upper,
                                                                                   x_rho_theta,
                                                                                   x_rho_w,
                                                                                   schur_element,
                                                                                   gravity,
                                                                                   theta,
                                                                                   dp_drho_theta,
                                                                                   metric_factor,
                                                                                   derivative_hat,
                                                                                   speed_of_sound,
                                                                                   inv_speed_of_sound,
                                                                                   inv_boundary_weight,
                                                                                   gamma_dt,
                                                                                   ::Val{NNODES},
                                                                                   column,
                                                                                   layer) where {
                                                                                                 NNODES
                                                                                                 }
    a12_1, inv_a22_1, b_lower, c_lower = first_node_entries(theta, dp_drho_theta,
                                                            metric_factor,
                                                            speed_of_sound,
                                                            inv_speed_of_sound,
                                                            inv_boundary_weight,
                                                            gamma_dt,
                                                            Val(NNODES), column, layer)
    theta_M = theta[column, NNODES, layer]
    dp_drho_theta_M = dp_drho_theta[column, NNODES, layer]
    zero_ = zero(eltype(x_rho_w))

    # Columns of `B` of `rho_w_M` (`r = 1`) and `rho_theta_M` (`r = 2`) of this element
    for i in 1:(NNODES - 1)
        derivative_iM = metric_factor[column, i, layer] * derivative_hat[i, NNODES]
        # `b_rho = A13[:, M]`, `b_rho_theta = A23[:, M]`, and `b_rho_w = 0`
        x_rho_theta[i, 1] = derivative_iM * theta_M
        x_rho_w[i, 1] = -gravity[column, i, layer] * gamma_dt * derivative_iM
        # `b_rho = 0`, `b_rho_theta = 0`, and `b_rho_w = A32[:, M]`
        x_rho_theta[i, 2] = zero_
        x_rho_w[i, 2] = derivative_iM * dp_drho_theta_M
    end
    eliminate_local_one_sided!(x_rho_theta, x_rho_w, schur_element, gravity, theta,
                               dp_drho_theta, metric_factor, derivative_hat, a12_1,
                               inv_a22_1, gamma_dt, Val(NNODES), column, layer)
    for r in 1:2
        minus_c_x = minus_c_x_one_sided(x_rho_theta, x_rho_w, theta, dp_drho_theta,
                                        metric_factor, derivative_hat, Val(NNODES),
                                        column, layer, r)
        for i in 1:3
            schur_diagonal[column, i, r + 1, layer] += minus_c_x[i]
        end
        # `C` of the layer below couples to the first node of this layer
        if layer > 1
            for i in 1:3
                schur_upper[column, i, r, layer - 1] -= c_lower[i, 1] * x_rho_w[1, r] +
                                                        c_lower[i, 2] *
                                                        x_rho_theta[1, r]
            end
        end
    end

    # Columns of `B` of `rho_w_M` (`r = 1`) and `rho_theta_M` (`r = 2`) of the layer below,
    # which are nonzero only in the first node
    if layer > 1
        gravity_1 = gravity[column, 1, layer]
        for r in 1:2
            for i in 1:(NNODES - 1)
                x_rho_theta[i, r] = zero_
                x_rho_w[i, r] = zero_
            end
            # `b_rho_theta` and `b_rho_w - A31 A11^{-1} b_rho`
            x_rho_theta[1, r] = b_lower[3, r]
            x_rho_w[1, r] = b_lower[2, r] - gravity_1 * gamma_dt * b_lower[1, r]
        end
        eliminate_local_one_sided!(x_rho_theta, x_rho_w, schur_element, gravity, theta,
                                   dp_drho_theta, metric_factor, derivative_hat, a12_1,
                                   inv_a22_1, gamma_dt, Val(NNODES), column, layer)
        for r in 1:2
            minus_c_x = minus_c_x_one_sided(x_rho_theta, x_rho_w, theta, dp_drho_theta,
                                            metric_factor, derivative_hat, Val(NNODES),
                                            column, layer, r)
            for i in 1:3
                schur_lower[column, i, r, layer - 1] += minus_c_x[i]
                # `C` of the layer below couples to the first node of this layer
                schur_diagonal[column, i, r + 1, layer - 1] -= c_lower[i, 1] *
                                                               x_rho_w[1, r] +
                                                               c_lower[i, 2] *
                                                               x_rho_theta[1, r]
            end
        end
    end
    return nothing
end

@kernel inbounds=true function schur_complement_boundary_one_sided_KAkernel!(schur_diagonal,
                                                                             schur_lower,
                                                                             schur_upper,
                                                                             @Const(schur_element),
                                                                             @Const(gravity),
                                                                             @Const(theta),
                                                                             @Const(dp_drho_theta),
                                                                             @Const(metric_factor),
                                                                             @Const(derivative_hat),
                                                                             speed_of_sound,
                                                                             inv_speed_of_sound,
                                                                             inv_boundary_weight,
                                                                             gamma_dt,
                                                                             ::Val{NLAYERS},
                                                                             ::Val{NNODES}) where {
                                                                                                   NLAYERS,
                                                                                                   NNODES
                                                                                                   }
    column = @index(Global)
    x_rho_theta = MMatrix{NNODES - 1, 2, eltype(schur_diagonal)}(undef)
    x_rho_w = MMatrix{NNODES - 1, 2, eltype(schur_diagonal)}(undef)

    for layer in 1:NLAYERS
        schur_complement_boundary_one_sided_per_element!(schur_diagonal, schur_lower,
                                                         schur_upper, x_rho_theta,
                                                         x_rho_w, schur_element,
                                                         gravity,
                                                         theta, dp_drho_theta,
                                                         metric_factor, derivative_hat,
                                                         speed_of_sound,
                                                         inv_speed_of_sound,
                                                         inv_boundary_weight, gamma_dt,
                                                         Val(NNODES), column, layer)
    end
end

function schur_complement_boundary!(backend::Backend,
                                    jacobian::JacobianVerticalOneSidedCache,
                                    gamma_dt)
    (; schur_diagonal, schur_lower, schur_upper, schur_element, gravity, theta,
    dp_drho_theta, metric_factor, derivative_hat, inv_boundary_weight, ncolumns) = jacobian
    speed_of_sound = jacobian.speed_of_sound
    inv_speed_of_sound = inv(speed_of_sound)
    nnodes_ = nnodes(jacobian)
    nlayers_ = nlayers(jacobian)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck check_axes_columns(jacobian)

    kernel! = schur_complement_boundary_one_sided_KAkernel!(backend)
    kernel!(schur_diagonal, schur_lower, schur_upper, schur_element, gravity, theta,
            dp_drho_theta, metric_factor, derivative_hat, speed_of_sound,
            inv_speed_of_sound, inv_boundary_weight, gamma_dt, Val(nlayers_),
            Val(nnodes_);
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Block LU factorization of the boundary Schur complement and the corresponding solve.
# The edge DOF `rho_M` (`u`) of a layer couples only to `rho_w_M` and `rho_theta_M` (`v`)
# of the same layer, since the first column of the blocks off the diagonal vanishes, and
# `D[3, 1] = 0`. Thus, eliminating `u` from the second row of each diagonal block with the
# multiplier `m = D[2, 1] / D[1, 1]` leaves a block-tridiagonal system with 2x2 blocks for
# `v`, which is factorized by a block LU without pivoting:
#   D_1' = D_1, L_i' = L_i D_{i-1}'^{-1}, D_i' = D_i - L_i' U_{i-1},
# where `D_i`, `L_i`, and `U_i` are the 2x2 blocks of `v` after the elimination of `u`.
# Afterward, `u` follows from the first row of each diagonal block.

# 2x2 blocks of `v` of the diagonal block `layer`, of the block below the diagonal in the
# block row `layer` (i.e., `schur_lower[:, :, :, layer - 1]`), and of the block above the
# diagonal in the block row `layer`, after the elimination of `u` with the multiplier `m`
Base.@propagate_inbounds function block_diagonal_one_sided(schur_diagonal, m, column,
                                                           layer)
    return SMatrix{2, 2}(schur_diagonal[column, 2, 2, layer] -
                         m * schur_diagonal[column, 1, 2, layer],
                         schur_diagonal[column, 3, 2, layer],
                         schur_diagonal[column, 2, 3, layer] -
                         m * schur_diagonal[column, 1, 3, layer],
                         schur_diagonal[column, 3, 3, layer])
end

Base.@propagate_inbounds function block_offdiagonal_one_sided(schur_offdiagonal, m,
                                                              column, index)
    return SMatrix{2, 2}(schur_offdiagonal[column, 2, 1, index] -
                         m * schur_offdiagonal[column, 1, 1, index],
                         schur_offdiagonal[column, 3, 1, index],
                         schur_offdiagonal[column, 2, 2, index] -
                         m * schur_offdiagonal[column, 1, 2, index],
                         schur_offdiagonal[column, 3, 2, index])
end

Base.@propagate_inbounds function store_block!(A, block, column, index)
    for j in 1:2, i in 1:2
        A[column, i, j, index] = block[i, j]
    end
    return nothing
end

Base.@propagate_inbounds function load_block(A, column, index)
    return SMatrix{2, 2}(A[column, 1, 1, index], A[column, 2, 1, index],
                         A[column, 1, 2, index], A[column, 2, 2, index])
end

Base.@propagate_inbounds function lu_block_per_column!(lu_multiplier, lu_inv_pivot,
                                                       lu_inv_diagonal, lu_lower,
                                                       lu_upper, schur_diagonal,
                                                       schur_lower, schur_upper,
                                                       ::Val{NLAYERS},
                                                       column) where {NLAYERS}
    for layer in 1:NLAYERS
        inv_pivot = inv(schur_diagonal[column, 1, 1, layer])
        lu_inv_pivot[column, layer] = inv_pivot
        lu_multiplier[column, layer] = schur_diagonal[column, 2, 1, layer] * inv_pivot
    end

    m = lu_multiplier[column, 1]
    inv_diagonal = inv(block_diagonal_one_sided(schur_diagonal, m, column, 1))
    store_block!(lu_inv_diagonal, inv_diagonal, column, 1)
    for layer in 2:NLAYERS
        # The block above the diagonal in the block row `layer - 1`
        upper = block_offdiagonal_one_sided(schur_upper,
                                            lu_multiplier[column, layer - 1],
                                            column, layer - 1)
        store_block!(lu_upper, upper, column, layer - 1)

        m = lu_multiplier[column, layer]
        lower = block_offdiagonal_one_sided(schur_lower, m, column, layer - 1) *
                inv_diagonal
        store_block!(lu_lower, lower, column, layer - 1)
        inv_diagonal = inv(block_diagonal_one_sided(schur_diagonal, m, column, layer) -
                           lower * upper)
        store_block!(lu_inv_diagonal, inv_diagonal, column, layer)
    end
    return nothing
end

@kernel inbounds=true function lu_block_KAkernel!(lu_multiplier, lu_inv_pivot,
                                                  lu_inv_diagonal, lu_lower, lu_upper,
                                                  @Const(schur_diagonal),
                                                  @Const(schur_lower),
                                                  @Const(schur_upper),
                                                  ::Val{NLAYERS}) where {NLAYERS}
    column = @index(Global)

    lu_block_per_column!(lu_multiplier, lu_inv_pivot, lu_inv_diagonal, lu_lower,
                         lu_upper, schur_diagonal, schur_lower, schur_upper,
                         Val(NLAYERS), column)
end

function lu_block!(backend::Backend, jacobian::JacobianVerticalOneSidedCache)
    (; lu_multiplier, lu_inv_pivot, lu_inv_diagonal, lu_lower, lu_upper,
    schur_diagonal, schur_lower, schur_upper, ncolumns) = jacobian
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck check_axes_columns(jacobian)
    kernel! = lu_block_KAkernel!(backend)
    kernel!(lu_multiplier, lu_inv_pivot, lu_inv_diagonal, lu_lower, lu_upper,
            schur_diagonal, schur_lower, schur_upper, Val(nlayers(jacobian));
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

# Solve with the block LU factors and the right-hand side `b = schur_rhs`, which holds the
# edge DOFs `(u, v)` of every layer in the rows `3 * (layer - 1) .+ (1:3)`
Base.@propagate_inbounds function ldiv_block_per_column!(b, lu_multiplier, lu_inv_pivot,
                                                         lu_inv_diagonal, lu_lower,
                                                         lu_upper, schur_diagonal,
                                                         schur_lower, schur_upper,
                                                         ::Val{NLAYERS},
                                                         column) where {NLAYERS}
    # Elimination of `u` and forward substitution
    v_previous = zero(SVector{2, eltype(b)})
    for layer in 1:NLAYERS
        offset = 3 * (layer - 1)
        b[column, offset + 2] -= lu_multiplier[column, layer] * b[column, offset + 1]
        v = SVector(b[column, offset + 2], b[column, offset + 3])
        if layer > 1
            v -= load_block(lu_lower, column, layer - 1) * v_previous
        end
        b[column, offset + 2] = v[1]
        b[column, offset + 3] = v[2]
        v_previous = v
    end

    # Backward substitution
    for layer in NLAYERS:-1:1
        offset = 3 * (layer - 1)
        v = SVector(b[column, offset + 2], b[column, offset + 3])
        if layer < NLAYERS
            v -= load_block(lu_upper, column, layer) *
                 SVector(b[column, offset + 5], b[column, offset + 6])
        end
        v = load_block(lu_inv_diagonal, column, layer) * v
        b[column, offset + 2] = v[1]
        b[column, offset + 3] = v[2]
    end

    # `u` from the first row of the diagonal blocks
    for layer in 1:NLAYERS
        offset = 3 * (layer - 1)
        u = b[column, offset + 1] -
            schur_diagonal[column, 1, 2, layer] * b[column, offset + 2] -
            schur_diagonal[column, 1, 3, layer] * b[column, offset + 3]
        if layer < NLAYERS
            u -= schur_upper[column, 1, 1, layer] * b[column, offset + 5] +
                 schur_upper[column, 1, 2, layer] * b[column, offset + 6]
        end
        if layer > 1
            u -= schur_lower[column, 1, 1, layer - 1] * b[column, offset - 1] +
                 schur_lower[column, 1, 2, layer - 1] * b[column, offset]
        end
        b[column, offset + 1] = lu_inv_pivot[column, layer] * u
    end
    return nothing
end

@kernel inbounds=true function ldiv_block_KAkernel!(schur_rhs, @Const(lu_multiplier),
                                                    @Const(lu_inv_pivot),
                                                    @Const(lu_inv_diagonal),
                                                    @Const(lu_lower), @Const(lu_upper),
                                                    @Const(schur_diagonal),
                                                    @Const(schur_lower),
                                                    @Const(schur_upper),
                                                    ::Val{NLAYERS}) where {NLAYERS}
    column = @index(Global)

    ldiv_block_per_column!(schur_rhs, lu_multiplier, lu_inv_pivot, lu_inv_diagonal,
                           lu_lower, lu_upper, schur_diagonal, schur_lower, schur_upper,
                           Val(NLAYERS), column)
end

function ldiv_block!(backend::Backend, jacobian::JacobianVerticalOneSidedCache)
    (; schur_rhs, lu_multiplier, lu_inv_pivot, lu_inv_diagonal, lu_lower, lu_upper,
    schur_diagonal, schur_lower, schur_upper, ncolumns) = jacobian
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck check_axes_columns(jacobian)
    kernel! = ldiv_block_KAkernel!(backend)
    kernel!(schur_rhs, lu_multiplier, lu_inv_pivot, lu_inv_diagonal, lu_lower, lu_upper,
            schur_diagonal, schur_lower, schur_upper, Val(nlayers(jacobian));
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Forward substitution: the right-hand side of the condensed system
#   `schur_rhs = b_D - C * A^{-1} * b_A`,
# where `A` is the block of the nodes `1..M-1` and `D` the block of the edge degrees of
# freedom

# Entries of `schur_rhs` of element `layer` of the column `column`, and the contribution
# `-C x` of this element to the edge DOFs of the layer below in `schur_rhs_lower`, which is
# added by `add_schur_rhs_lower_KAkernel!`. `x_rho_theta` and `x_rho_w` hold the components
# of `A^{-1} b_A`, see `eliminate_local_one_sided!`.
Base.@propagate_inbounds function ldiv_interior_forward_one_sided_per_element!(schur_rhs,
                                                                               schur_rhs_lower,
                                                                               x_rho_theta,
                                                                               x_rho_w,
                                                                               b,
                                                                               schur_element,
                                                                               gravity,
                                                                               theta,
                                                                               dp_drho_theta,
                                                                               metric_factor,
                                                                               derivative_hat,
                                                                               speed_of_sound,
                                                                               inv_speed_of_sound,
                                                                               inv_boundary_weight,
                                                                               gamma_dt,
                                                                               ::Val{NNODES},
                                                                               column,
                                                                               layer) where {
                                                                                             NNODES
                                                                                             }
    a12_1, inv_a22_1, _, c_lower = first_node_entries(theta, dp_drho_theta,
                                                      metric_factor, speed_of_sound,
                                                      inv_speed_of_sound,
                                                      inv_boundary_weight, gamma_dt,
                                                      Val(NNODES), column, layer)
    for i in 1:(NNODES - 1)
        x_rho_theta[i] = b[column, 3, i, layer]
        x_rho_w[i] = b[column, 2, i, layer] -
                     gravity[column, i, layer] * gamma_dt * b[column, 1, i, layer]
    end
    eliminate_local_one_sided!(x_rho_theta, x_rho_w, schur_element, gravity, theta,
                               dp_drho_theta, metric_factor, derivative_hat, a12_1,
                               inv_a22_1, gamma_dt, Val(NNODES), column, layer)
    minus_c_x = minus_c_x_one_sided(x_rho_theta, x_rho_w, theta, dp_drho_theta,
                                    metric_factor, derivative_hat, Val(NNODES), column,
                                    layer, 1)
    offset = 3 * (layer - 1)
    for i in 1:3
        schur_rhs[column, offset + i] = b[column, i, NNODES, layer] + minus_c_x[i]
    end
    if layer > 1
        for i in 1:3
            schur_rhs_lower[column, i, layer] = -(c_lower[i, 1] * x_rho_w[1] +
                                                  c_lower[i, 2] * x_rho_theta[1])
        end
    end
    return nothing
end

@kernel inbounds=true function ldiv_interior_forward_one_sided_KAkernel!(schur_rhs,
                                                                         schur_rhs_lower,
                                                                         @Const(b),
                                                                         @Const(schur_element),
                                                                         @Const(gravity),
                                                                         @Const(theta),
                                                                         @Const(dp_drho_theta),
                                                                         @Const(metric_factor),
                                                                         @Const(derivative_hat),
                                                                         speed_of_sound,
                                                                         inv_speed_of_sound,
                                                                         inv_boundary_weight,
                                                                         gamma_dt,
                                                                         ::Val{NNODES}) where {
                                                                                               NNODES
                                                                                               }
    column, layer = @index(Global, NTuple)
    x_rho_theta = MVector{NNODES - 1, eltype(b)}(undef)
    x_rho_w = MVector{NNODES - 1, eltype(b)}(undef)

    ldiv_interior_forward_one_sided_per_element!(schur_rhs, schur_rhs_lower,
                                                 x_rho_theta,
                                                 x_rho_w, b, schur_element, gravity,
                                                 theta, dp_drho_theta, metric_factor,
                                                 derivative_hat, speed_of_sound,
                                                 inv_speed_of_sound,
                                                 inv_boundary_weight,
                                                 gamma_dt, Val(NNODES), column, layer)
end

# Add the contributions of the layer above to the edge DOFs of `layer`
@kernel inbounds=true function add_schur_rhs_lower_KAkernel!(schur_rhs,
                                                             @Const(schur_rhs_lower))
    column, layer = @index(Global, NTuple)
    offset = 3 * (layer - 1)
    for i in 1:3
        schur_rhs[column, offset + i] += schur_rhs_lower[column, i, layer + 1]
    end
end

function ldiv_interior_forward!(backend::Backend,
                                jacobian::JacobianVerticalOneSidedCache, b,
                                gamma_dt)
    (; schur_rhs, schur_rhs_lower, schur_element, gravity, theta, dp_drho_theta,
    metric_factor, derivative_hat, inv_boundary_weight, ncolumns) = jacobian
    speed_of_sound = jacobian.speed_of_sound
    inv_speed_of_sound = inv(speed_of_sound)
    nnodes_ = nnodes(jacobian)
    nlayers_ = nlayers(jacobian)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        check_axes_columns(jacobian)
        Trixi.check_axes(b, axes(jacobian.u_vertical))
    end

    kernel! = ldiv_interior_forward_one_sided_KAkernel!(backend)
    kernel!(schur_rhs, schur_rhs_lower, b, schur_element, gravity, theta, dp_drho_theta,
            metric_factor, derivative_hat, speed_of_sound, inv_speed_of_sound,
            inv_boundary_weight, gamma_dt, Val(nnodes_);
            ndrange = (ncolumns, nlayers_), layer_kernel_workgroupsize()...)
    if nlayers_ > 1
        kernel! = add_schur_rhs_lower_KAkernel!(backend)
        kernel!(schur_rhs, schur_rhs_lower;
                ndrange = (ncolumns, nlayers_ - 1), layer_kernel_workgroupsize()...)
    end
    return nothing
end

###############################################################################
# Backward substitution: the solution of the edge degrees of freedom `schur_rhs` and of
# the degrees of freedom in the nodes `1..M-1`

# Solution of element `layer` of the column `column`; `x_rho`, `x_rho_theta`, and `x_rho_w`
# hold the components of the three variables in the nodes `1..M-1`
Base.@propagate_inbounds function ldiv_interior_backward_one_sided_per_element!(b,
                                                                                x_rho,
                                                                                x_rho_theta,
                                                                                x_rho_w,
                                                                                schur_rhs,
                                                                                schur_element,
                                                                                gravity,
                                                                                theta,
                                                                                dp_drho_theta,
                                                                                metric_factor,
                                                                                derivative_hat,
                                                                                speed_of_sound,
                                                                                inv_speed_of_sound,
                                                                                inv_boundary_weight,
                                                                                gamma_dt,
                                                                                ::Val{NNODES},
                                                                                column,
                                                                                layer) where {
                                                                                              NNODES
                                                                                              }
    a12_1, inv_a22_1, b_lower, _ = first_node_entries(theta, dp_drho_theta,
                                                      metric_factor, speed_of_sound,
                                                      inv_speed_of_sound,
                                                      inv_boundary_weight, gamma_dt,
                                                      Val(NNODES), column, layer)
    offset = 3 * (layer - 1)
    rho_w_M = schur_rhs[column, offset + 2]
    rho_theta_M = schur_rhs[column, offset + 3]
    for v in 1:3
        b[column, v, NNODES, layer] = schur_rhs[column, offset + v]
    end

    # `b_A - B x_D` with the columns of `rho_w_M` and `rho_theta_M` of this element
    theta_M = theta[column, NNODES, layer]
    dp_drho_theta_M = dp_drho_theta[column, NNODES, layer]
    for i in 1:(NNODES - 1)
        derivative_iM = metric_factor[column, i, layer] * derivative_hat[i, NNODES]
        x_rho[i] = b[column, 1, i, layer] - derivative_iM * rho_w_M
        x_rho_w[i] = b[column, 2, i, layer] -
                     derivative_iM * dp_drho_theta_M * rho_theta_M
        x_rho_theta[i] = b[column, 3, i, layer] - derivative_iM * theta_M * rho_w_M
    end
    # and with the columns of `rho_w_M` and `rho_theta_M` of the layer below
    if layer > 1
        rho_w_ll = schur_rhs[column, offset - 1]
        rho_theta_ll = schur_rhs[column, offset]
        x_rho[1] -= b_lower[1, 1] * rho_w_ll + b_lower[1, 2] * rho_theta_ll
        x_rho_w[1] -= b_lower[2, 1] * rho_w_ll + b_lower[2, 2] * rho_theta_ll
        x_rho_theta[1] -= b_lower[3, 1] * rho_w_ll + b_lower[3, 2] * rho_theta_ll
    end

    # `x_A = A^{-1} (b_A - B x_D)`
    for i in 1:(NNODES - 1)
        x_rho_w[i] -= gravity[column, i, layer] * gamma_dt * x_rho[i]
    end
    eliminate_local_one_sided!(x_rho_theta, x_rho_w, schur_element, gravity, theta,
                               dp_drho_theta, metric_factor, derivative_hat, a12_1,
                               inv_a22_1, gamma_dt, Val(NNODES), column, layer)
    # x_rho = A11^{-1} (b_rho - A12 x_rho_theta - A13 x_rho_w)
    factor_11_a13 = first(factors_11(layer))
    x_rho[1] -= a12_1 * x_rho_theta[1]
    for j in 1:(NNODES - 1)
        x_rho_w_j = x_rho_w[j]
        for i in 1:(NNODES - 1)
            x_rho[i] -= derivative_entry(derivative_hat,
                                         metric_factor[column, i, layer],
                                         factor_11_a13, i, j) * x_rho_w_j
        end
    end

    for i in 1:(NNODES - 1)
        b[column, 1, i, layer] = gamma_dt * x_rho[i]
        b[column, 2, i, layer] = x_rho_w[i]
        b[column, 3, i, layer] = x_rho_theta[i]
    end
    return nothing
end

@kernel inbounds=true function ldiv_interior_backward_one_sided_KAkernel!(b,
                                                                          @Const(schur_rhs),
                                                                          @Const(schur_element),
                                                                          @Const(gravity),
                                                                          @Const(theta),
                                                                          @Const(dp_drho_theta),
                                                                          @Const(metric_factor),
                                                                          @Const(derivative_hat),
                                                                          speed_of_sound,
                                                                          inv_speed_of_sound,
                                                                          inv_boundary_weight,
                                                                          gamma_dt,
                                                                          ::Val{NNODES}) where {
                                                                                                NNODES
                                                                                                }
    column, layer = @index(Global, NTuple)
    x_rho = MVector{NNODES - 1, eltype(b)}(undef)
    x_rho_theta = MVector{NNODES - 1, eltype(b)}(undef)
    x_rho_w = MVector{NNODES - 1, eltype(b)}(undef)

    ldiv_interior_backward_one_sided_per_element!(b, x_rho, x_rho_theta, x_rho_w,
                                                  schur_rhs, schur_element, gravity,
                                                  theta, dp_drho_theta, metric_factor,
                                                  derivative_hat, speed_of_sound,
                                                  inv_speed_of_sound,
                                                  inv_boundary_weight,
                                                  gamma_dt, Val(NNODES), column, layer)
end

function ldiv_interior_backward!(backend::Backend,
                                 jacobian::JacobianVerticalOneSidedCache,
                                 b,
                                 gamma_dt)
    (; schur_rhs, schur_element, gravity, theta, dp_drho_theta, metric_factor,
    derivative_hat, inv_boundary_weight, ncolumns) = jacobian
    speed_of_sound = jacobian.speed_of_sound
    inv_speed_of_sound = inv(speed_of_sound)
    nnodes_ = nnodes(jacobian)
    nlayers_ = nlayers(jacobian)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        check_axes_columns(jacobian)
        Trixi.check_axes(b, axes(jacobian.u_vertical))
    end

    kernel! = ldiv_interior_backward_one_sided_KAkernel!(backend)
    kernel!(b, schur_rhs, schur_element, gravity, theta, dp_drho_theta, metric_factor,
            derivative_hat, speed_of_sound, inv_speed_of_sound, inv_boundary_weight,
            gamma_dt, Val(nnodes_);
            ndrange = (ncolumns, nlayers_), layer_kernel_workgroupsize()...)
    return nothing
end
end # @muladd
