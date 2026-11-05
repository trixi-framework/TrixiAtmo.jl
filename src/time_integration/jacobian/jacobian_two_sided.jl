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
arrays are created by [`create_jacobian_cache`](@ref) when the time integrator is
initialized, for the semidiscretization of the ODE problem and with the array type of the
solution, so that it runs on the same backend (CPU or GPU) as the semidiscretization. The
Jacobian decouples into independent linear systems for each vertical column of nodes. The
dissipation of the vertical acoustic waves at the element interfaces and at the walls is the
one of the LMARS flux with the constant speed of sound `c`, also if the surface flux of the
semidiscretization uses a local speed of sound (e.g., [`flux_lmars_es`](@ref)), since the
Jacobian only needs to approximate the right-hand side.

### Permutation
The unknowns of a column are the density `rho`, the vertical momentum `rho_w`, i.e., the
component of the momentum in the unit contravariant direction of the vertical reference
coordinate, and the potential temperature density `rho_theta`. The horizontal components
of the momentum are not part of the Jacobian. For a column of `nz` elements of polynomial
degree `polydeg` (`M = polydeg + 1` LGL nodes), the unknowns are grouped, by variable, into
a small *local* block per element and a small *global* (edge) block shared by the whole
column:

* all `rho` DOFs of element 1 (nodes `1..M`), then
* the `rho_theta` DOFs `2..M-1` of element 1 (interior only), then
* the `rho_w` DOFs `2..M-1` of element 1 (interior only), then
* ... (repeated, element by element) ...
* finally, 4 edge DOFs per element, ordered element by element:
  `rho_theta_1, rho_w_1, rho_w_M, rho_theta_M`.

Density is *never* condensed: it keeps all `M` nodes, including both edges,
in the local block. Only `rho_theta` and `rho_w` are condensed, and only at
*both* edge nodes of each element -- that is what "two-sided" refers to.

### Sparsity pattern
Example for `nz = 4`, `polydeg = 5` (`M = 6`), in the permuted ordering above:
```
⠑⢄⠀⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀
⠀⠀⠑⢄⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡗⠒⠀⠀⠀⠀⠀
⠠⡀⠀⣤⣵⢟⠛⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⢣⠀⠀⠀⠀⠀⠀
⠀⠈⠂⠛⠛⠀⠑⢄⠀⠀⠀⠀⣤⣤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠃⠼⢤⡄⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣇⣀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠈⠢⡀⣿⣿⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⢸⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢹⡇⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡗⠒⠀ = ⎡A B⎤
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⡀⠀⣤⣵⢟⠛⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⢣⠀⠀   ⎣C D⎦
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠂⠛⠛⠀⠑⢄⠀⠀⠀⠀⣤⣤⠀⠀⠀⠀⠃⠼⢤⡄
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⢸⡇
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⣿⣿⠀⠀⠀⠀⠀⠀⢸⡇
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⡀⣿⣿⠑⢄⠀⠀⠀⠀⠀⠀⡇⢸
⠂⠀⠠⠶⠶⣉⣉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢟⣵⣤⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠂⠀⠠⠶⠶⣉⣉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⢟⣵⣤⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠠⠶⠶⣉⣉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⢟⣵⣤⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠠⠶⠶⣉⣉⠀⠀⠀⠀⠀⠛⢟⣵
```

`A` (block-diagonal, one block per element) holds the diagonals `I / gamma_dt` of `rho`
and of the interior `rho_theta`/`rho_w`, the dense coupling (LGL derivative operator) of
`rho` and the interior `rho_theta` to the interior `rho_w` (mass and potential temperature
fluxes) and of the interior `rho_w` to the interior `rho_theta` (pressure gradient), and
the pointwise coupling of the interior `rho_w` to `rho` (gravity):
```
⠑⢄⠀⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠑⢄⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠠⡀⠀⣤⣵⢟⠛⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠈⠂⠛⠛⠀⠑⢄⠀⠀⠀⠀⣤⣤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠈⠢⡀⣿⣿⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⡀⠀⣤⣵⢟⠛⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠂⠛⠛⠀⠑⢄⠀⠀⠀⠀⣤⣤
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⣿⣿
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⣿⣿
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⡀⣿⣿⠑⢄
```

`B` connects the local DOFs to both edges of `rho_theta`/`rho_w` of the same
element (dense LGL coupling) and, through the LMARS mass flux, the edge `rho`
to the edge `rho_theta` of the same element and to the edge DOFs of the neighboring
elements:
```
⢸⡇⠀⠀⠀⠀⠀⠀
⢸⡗⠒⠀⠀⠀⠀⠀
⡜⢣⠀⠀⠀⠀⠀⠀
⠃⠼⢤⡄⠀⠀⠀⠀
⠀⠀⢸⣇⣀⠀⠀⠀
⠀⠀⢸⡇⠀⠀⠀⠀
⠀⠀⡇⢸⠀⠀⠀⠀
⠀⠀⠀⠉⢹⡇⠀⠀
⠀⠀⠀⠀⢸⡗⠒⠀
⠀⠀⠀⠀⡜⢣⠀⠀
⠀⠀⠀⠀⠃⠼⢤⡄
⠀⠀⠀⠀⠀⠀⢸⡇
⠀⠀⠀⠀⠀⠀⢸⡇
⠀⠀⠀⠀⠀⠀⡇⢸
```

`C` connects the edge `rho_theta`/`rho_w` to the interior `rho_w`/`rho_theta` of
the same element (dense LGL coupling) and, pointwise, the edge `rho_w` to the
local `rho` (gravity):
```
⠂⠀⠠⠶⠶⣉⣉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠂⠀⠠⠶⠶⣉⣉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠠⠶⠶⣉⣉⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠠⠶⠶⣉⣉
```

`D` is block-tridiagonal before elimination: the 4 edge DOFs of each element couple to
each other (dense LGL coupling and the LMARS flux) and to those of the neighboring
elements through the single-valued LMARS flux:
```
⢟⣵⣤⠀⠀⠀⠀⠀
⠀⠛⢟⣵⣤⠀⠀⠀
⠀⠀⠀⠛⢟⣵⣤⠀
⠀⠀⠀⠀⠀⠛⢟⣵
```

### Algorithm
The matrix `I / gamma_dt - J` is assembled and factorized by `update_jacobian!`, i.e.,
whenever the Jacobian is updated (see the keyword argument `update_interval` of
`Trixi.solve`), and the linear systems of the stages are solved by `solve_jacobian!`:
1. **Local Schur complement** (per element, in `update_jacobian!`): eliminate the
   diagonal blocks of `rho` and the interior `rho_theta` analytically, leaving a dense
   `(M-2)x(M-2)` complement `S_e` on the interior `rho_w`, which is LU-factorized
   without pivoting (stored in `schur_element`).
2. **Global Schur complement** (per column, in `update_jacobian!`): eliminating the
   local DOFs of every element with `S_e` leaves a `4*nz x 4*nz` banded system
   `S = D - C A^{-1} B` (3 sub- and 3 super-diagonals) on the edge DOFs, which is
   factorized with a banded LU without pivoting (stored in `schur_band`).
3. **Solve** (per stage, in `solve_jacobian!`): copy the right-hand side to the columns,
   eliminate the local DOFs with `S_e` (right-hand side `b_D - C A^{-1} b_A` of the
   condensed system), solve the banded system with the LU factors of `S` for the edge
   DOFs, back-substitute into the local DOFs, and copy the solution back to Trixi.jl's
   layout. The horizontal components of the momentum are only multiplied by `gamma_dt`,
   i.e., solved with `J = 0`.

All kernels use KernelAbstractions.jl with the backend of the arrays of the
semidiscretization, i.e., the CPU backend of KernelAbstractions.jl for `Array`s
(independently of the threading backend of Trixi.jl) or a GPU backend. Columns are
independent, so the sequential recurrences (fill, Schur assembly, banded factorization
and solve) use one thread per column, the copies between Trixi.jl's layout and the columns
use one thread per node, and the independent work per element (forward elimination,
back-substitution) uses one thread per `(column, element)`.

!!! note
    The mesh must not be refined, since the elements of a column must be consecutive
    with the lowest layer first. For parallel meshes, the columns must be kept on one rank
    (`keep_columns_together = true` in [`P4estMeshCubedSphereTopography`](@ref)).
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
    dp_drho_theta::ArrayRealT3
    # LU factors of the local Schur complements `S_e` of the elements in the interior
    # nodes, `[column, node, node, layer]`
    schur_element::ArrayRealT4
    # Global Schur complement `S`, a band matrix with 7 diagonals laid out
    # `[column, diagonal, row]`, and the right-hand side and solution of the condensed
    # system, `[column, row]`
    schur_band::ArrayRealT3
    schur_rhs::MatrixRealT
end

# Number of layers of elements of a column and number of nodes per coordinate direction of
# the elements (as `nnodes(solver)`), stored as type parameters such that they are
# compile-time constants in the kernels
@inline function nlayers(::VerticalJacobianTwoSidedCache{NLAYERS}) where {NLAYERS}
    return NLAYERS
end

@inline function Trixi.nnodes(::VerticalJacobianTwoSidedCache{<:Any, NNODES}) where {NNODES}
    return NNODES
end

# Check whether the arrays of the columns in `jacobian` have the axes we assume they must
# have in the kernels of the vertical Jacobian.
function check_axes_columns(jacobian::VerticalJacobianTwoSidedCache)
    (; ncolumns, metric_factor, gravity, derivative_hat, u_vertical, theta,
    dp_drho_theta, schur_element, schur_band, schur_rhs) = jacobian
    columns = Base.OneTo(ncolumns)
    nodes = Base.OneTo(nnodes(jacobian))
    interior_nodes = Base.OneTo(nnodes(jacobian) - 2)
    layers = Base.OneTo(nlayers(jacobian))
    Trixi.check_axes(metric_factor, (columns, nodes, layers))
    Trixi.check_axes(gravity, (columns, nodes, layers))
    Trixi.check_axes(theta, (columns, nodes, layers))
    Trixi.check_axes(dp_drho_theta, (columns, nodes, layers))
    Trixi.check_axes(derivative_hat, (nodes, nodes))
    Trixi.check_axes(u_vertical, (columns, Base.OneTo(3), nodes, layers))
    Trixi.check_axes(schur_element, (columns, interior_nodes, interior_nodes, layers))
    # Band matrix with 3 lower and 3 upper diagonals and 4 boundary unknowns per layer
    Trixi.check_axes(schur_band,
                     (columns, Base.OneTo(7), Base.OneTo(4 * length(layers))))
    Trixi.check_axes(schur_rhs, (columns, Base.OneTo(4 * length(layers))))
    return nothing
end

function create_jacobian_cache(jacobian::VerticalJacobianTwoSided, semi, u_ode)
    (; mesh, equations, solver, cache) = semi
    if !(equations isa CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
        throw(ArgumentError("`VerticalJacobianTwoSided` requires the " *
                            "`CompressibleEulerPotentialTemperatureEquationsWithGravity3D`"))
    end
    nlayers_ = nlayers(mesh)
    nnodes_ = nnodes(solver)
    nelements_ = nelements(solver, cache)
    # Compare the global numbers of elements and trees, since `nelements` is local to the
    # MPI rank for parallel meshes
    if Trixi.ncellsglobal(mesh) != Trixi.ntrees(mesh)
        throw(ArgumentError("`VerticalJacobianTwoSided` requires a mesh " *
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
    vertical_frame = Adapt.adapt(storage, compute_vertical_frame(semi))
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
    schur_element = KernelAbstractions.zeros(backend, RealT, ncolumns, nnodes_ - 2,
                                             nnodes_ - 2,
                                             nlayers_)
    schur_band = KernelAbstractions.zeros(backend, RealT, ncolumns, 7, 4 * nlayers_)
    schur_rhs = KernelAbstractions.zeros(backend, RealT, ncolumns, 4 * nlayers_)

    return VerticalJacobianTwoSidedCache{nlayers_, nnodes_, RealT, typeof(backend),
                                         typeof(schur_rhs), typeof(gravity),
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
                                                                 theta, dp_drho_theta,
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
# factors `schur_element` of the local Schur complements, and the direct contributions to
# the boundary Schur complement `schur_band`. The density entries of `B` are computed where
# they are needed.

# Fill of element `layer` of the column `column`; the layers of a column are filled in
# order by the same thread
Base.@propagate_inbounds function fill_jacobian_per_element!(schur_element, schur_band,
                                                             theta, dp_drho_theta,
                                                             theta_element,
                                                             dp_drho_theta_element,
                                                             u_vertical, metric_factor,
                                                             gravity, derivative_hat,
                                                             speed_of_sound,
                                                             inv_speed_of_sound,
                                                             inv_boundary_weight,
                                                             dp_factor, dp_exponent,
                                                             inv_gamma_dt, gamma_dt,
                                                             nlayers_, ::Val{NNODES},
                                                             column,
                                                             layer) where {NNODES}
    RealT = eltype(u_vertical)
    zero_ = zero(RealT)
    # Index of the first of the 4 boundary unknowns of this layer in the condensed system
    index_boundary = (layer - 1) * 4 + 1
    for i in 1:NNODES
        theta_element[i] = u_vertical[column, 3, i, layer] /
                           u_vertical[column, 1, i, layer]
        dp_drho_theta_element[i] = dp_factor *
                                   u_vertical[column, 3, i, layer]^dp_exponent
        theta[column, i, layer] = theta_element[i]
        dp_drho_theta[column, i, layer] = dp_drho_theta_element[i]
    end

    # Local Schur complement; the derivative blocks `A13`, `A23`, and `A32` are computed
    # from `derivative_hat`, `theta_element`, `dp_drho_theta_element`, and `metric_factor`
    for i in 1:(NNODES - 2)
        gravity_i = gravity[column, i + 1, layer]
        metric_factor_i = metric_factor[column, i + 1, layer]
        for j in 1:(NNODES - 2)
            # Entry `(A32 * A23)[i, j]`; the factors `metric_factor_i` of `A32` and
            # `theta_element[j + 1]` of `A23` do not depend on `k`
            correction = zero_
            for k in 1:(NNODES - 2)
                correction += derivative_hat[i + 1, k + 1] *
                              dp_drho_theta_element[k + 1] *
                              metric_factor[column, k + 1, layer] *
                              derivative_hat[k + 1, j + 1]
            end
            correction *= metric_factor_i * theta_element[j + 1]
            a13 = metric_factor_i * derivative_hat[i + 1, j + 1]
            schur_element[column, i, j, layer] = (i == j ? inv_gamma_dt : zero_) -
                                                 gamma_dt *
                                                 (gravity_i * a13 + correction)
        end
    end
    lu_dense!(schur_element, Val(NNODES - 2), column, layer)

    if layer > 1
        # States of the lower interface, with the last node of the layer below as `ll`
        # (stored by the same thread when filling that layer) and the first node of this
        # layer as `rr`
        theta_ll = theta[column, NNODES, layer - 1]
        dp_drho_theta_ll = dp_drho_theta[column, NNODES, layer - 1]
        theta_rr = theta_element[1]
        dp_drho_theta_rr = dp_drho_theta_element[1]
        # Interface potential temperature of the dissipative flux of `rho_theta`
        theta_interface = 0.5f0 * (theta_ll + theta_rr)
        # Factors of the rows of the `ll` node (rows `index_boundary - 2` and
        # `index_boundary - 1`) and of the `rr` node (rows `index_boundary` and
        # `index_boundary + 1`), as in `schur_complement_boundary_per_element!`
        factor_ll = 0.5f0 * metric_factor[column, NNODES, layer - 1] *
                    inv_boundary_weight
        factor_rr = 0.5f0 * metric_factor[column, 1, layer] * inv_boundary_weight
        # Contributions of the LMARS flux at the interface. Each comment gives the entry
        # `(i, j)` (row, column) of the condensed system, stored at
        # `schur_band[column, i - j + 4, j]`.
        # (index_boundary - 1, index_boundary - 1)
        schur_band[column, 4, index_boundary - 1] += theta_interface *
                                                     dp_drho_theta_ll *
                                                     inv_speed_of_sound * factor_ll
        # (index_boundary, index_boundary - 1)
        schur_band[column, 5, index_boundary - 1] += -theta_interface *
                                                     dp_drho_theta_ll *
                                                     inv_speed_of_sound * factor_rr
        # (index_boundary - 1, index_boundary)
        schur_band[column, 3, index_boundary] += -theta_interface * dp_drho_theta_rr *
                                                 inv_speed_of_sound * factor_ll
        # (index_boundary, index_boundary)
        schur_band[column, 4, index_boundary] += theta_interface * dp_drho_theta_rr *
                                                 inv_speed_of_sound * factor_rr
        # (index_boundary - 2, index_boundary - 2)
        schur_band[column, 4, index_boundary - 2] += speed_of_sound * factor_ll
        # (index_boundary + 1, index_boundary - 2)
        schur_band[column, 7, index_boundary - 2] += -speed_of_sound * factor_rr
        # (index_boundary - 2, index_boundary + 1)
        schur_band[column, 1, index_boundary + 1] += -speed_of_sound * factor_ll
        # (index_boundary + 1, index_boundary + 1)
        schur_band[column, 4, index_boundary + 1] += speed_of_sound * factor_rr
        # (index_boundary, index_boundary - 2)
        schur_band[column, 6, index_boundary - 2] += -theta_ll * factor_rr
        # (index_boundary + 1, index_boundary - 1)
        schur_band[column, 6, index_boundary - 1] += -dp_drho_theta_ll * factor_rr
        # (index_boundary - 1, index_boundary + 1)
        schur_band[column, 2, index_boundary + 1] += theta_rr * factor_ll
        # (index_boundary - 2, index_boundary)
        schur_band[column, 2, index_boundary] += dp_drho_theta_rr * factor_ll
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
                                                 dp_drho_theta_element[1]
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
                                                 dp_drho_theta_element[NNODES]
    schur_band[column, 6, index_boundary + 1] += metric_factor[column,
                                                               NNODES,
                                                               layer] *
                                                 derivative_hat[NNODES, 1] *
                                                 theta_element[1]
    schur_band[column, 6, index_boundary] += metric_factor[column, NNODES,
                                                           layer] *
                                             derivative_hat[NNODES, 1] *
                                             dp_drho_theta_element[1]
    if layer == nlayers_
        schur_band[column, 5, index_boundary + 2] += metric_factor[column,
                                                                   NNODES,
                                                                   layer] *
                                                     derivative_hat[NNODES, NNODES] *
                                                     theta_element[NNODES]
        schur_band[column, 3, index_boundary + 3] += -metric_factor[column,
                                                                    NNODES,
                                                                    layer] *
                                                     derivative_hat[NNODES, NNODES] *
                                                     dp_drho_theta_element[NNODES]
        schur_band[column, 4, index_boundary + 2] += metric_factor[column,
                                                                   NNODES,
                                                                   layer] *
                                                     speed_of_sound *
                                                     inv_boundary_weight
    end
    return nothing
end

@kernel inbounds=true function fill_jacobian_KAkernel!(schur_element, schur_band, theta,
                                                       dp_drho_theta,
                                                       @Const(u_vertical),
                                                       @Const(metric_factor),
                                                       @Const(gravity),
                                                       @Const(derivative_hat),
                                                       speed_of_sound,
                                                       inv_speed_of_sound,
                                                       inv_boundary_weight, dp_factor,
                                                       dp_exponent, inv_gamma_dt,
                                                       gamma_dt, nlayers_,
                                                       ::Val{NNODES}) where {NNODES}
    column = @index(Global)
    theta_element = MVector{NNODES, eltype(u_vertical)}(undef)
    dp_drho_theta_element = MVector{NNODES, eltype(u_vertical)}(undef)

    for layer in 1:nlayers_
        fill_jacobian_per_element!(schur_element, schur_band, theta, dp_drho_theta,
                                   theta_element, dp_drho_theta_element, u_vertical,
                                   metric_factor, gravity, derivative_hat,
                                   speed_of_sound, inv_speed_of_sound,
                                   inv_boundary_weight, dp_factor, dp_exponent,
                                   inv_gamma_dt, gamma_dt, nlayers_, Val(NNODES),
                                   column, layer)
    end
end

function fill_jacobian!(backend::Backend, jacobian::VerticalJacobianTwoSidedCache,
                        equations,
                        gamma_dt)
    (; schur_element, schur_band, gravity, theta, dp_drho_theta, ncolumns, metric_factor,
    derivative_hat, inv_boundary_weight, u_vertical) = jacobian
    speed_of_sound = jacobian.speed_of_sound
    nnodes_ = nnodes(jacobian)
    nlayers_ = nlayers(jacobian)
    dp_factor = equations.gamma * equations.K
    dp_exponent = equations.gamma - 1
    inv_speed_of_sound = inv(speed_of_sound)
    inv_gamma_dt = inv(gamma_dt)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck check_axes_columns(jacobian)
    # The fill and the boundary Schur complement accumulate into `schur_band`
    fill!(schur_band, zero(eltype(schur_band)))
    @views schur_band[:, 4, :] .= inv_gamma_dt

    kernel! = fill_jacobian_KAkernel!(backend)
    kernel!(schur_element, schur_band, theta, dp_drho_theta, u_vertical, metric_factor,
            gravity, derivative_hat, speed_of_sound, inv_speed_of_sound,
            inv_boundary_weight, dp_factor, dp_exponent, inv_gamma_dt, gamma_dt,
            nlayers_, Val(nnodes_);
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Boundary Schur complement: the contributions of the interior nodes of every element,
# eliminated with the local Schur complements `schur_element`, to the band matrix
# `schur_band`

# Contribution `-C A^{-1} B` of element `layer` of the column `column`; the layers of a
# column are processed in order by the same thread. Each of the 4 columns of `B`, i.e., of
# the boundary unknowns `rho_theta_1`, `rho_w_1`, `rho_w_M`, and `rho_theta_M` of the
# element (`index_local = 1, ..., 4`), is set up here and eliminated by
# `schur_complement_boundary_unknown!`. The arrays `x_rho`, `x_rho_theta`, and `x_rho_w`
# hold the components of the column `b` of `B` of the three variables.
Base.@propagate_inbounds function schur_complement_boundary_per_element!(schur_band,
                                                                         x_rho,
                                                                         x_rho_theta,
                                                                         x_rho_w,
                                                                         schur_element,
                                                                         gravity, theta,
                                                                         dp_drho_theta,
                                                                         metric_factor,
                                                                         derivative_hat,
                                                                         gamma_dt,
                                                                         inv_speed_of_sound,
                                                                         inv_boundary_weight,
                                                                         nlayers_,
                                                                         ::Val{NNODES},
                                                                         column,
                                                                         layer) where {NNODES}
    RealT = eltype(schur_band)
    zero_ = zero(RealT)
    # Offset of the 4 boundary unknowns of this layer in the condensed system
    offset_boundary = (layer - 1) * 4

    for index_local in 1:4
        # Edge node of the boundary unknown and whether it lies at an interface between
        # elements (the lower one for the first node, the upper one for the last node)
        # instead of at a wall
        if index_local <= 2
            node = 1
            is_interface = layer > 1
        else
            node = NNODES
            is_interface = layer < nlayers_
        end
        metric_factor_node = metric_factor[column, node, layer]

        # Components of `b` in `x_rho` and `x_rho_theta`, and the right-hand side
        #   b_rho_w - A31 A11^{-1} b_rho - A32 A22^{-1} b_rho_theta
        # of `S_e x_rho_w` in `x_rho_w` (see `schur_complement_boundary_unknown!`)
        if index_local == 1 || index_local == 4
            # Column of `rho_theta`: `b_rho_w = A32[:, node]`, `b_rho_theta = 0`, and the
            # density entry of `B` in the edge node from the LMARS mass flux at the
            # interface, which vanishes at the wall. Since `b_rho` vanishes in the interior
            # nodes, the right-hand side of `S_e x_rho_w` is `b_rho_w`, and only the edge
            # nodes of `x_rho` are used.
            x_rho[1] = zero_
            x_rho[NNODES] = zero_
            if is_interface
                factor = 0.5f0 * metric_factor_node * inv_boundary_weight
                x_rho[node] = dp_drho_theta[column, node, layer] * inv_speed_of_sound *
                              factor
            end
            for i in 1:(NNODES - 2)
                x_rho_theta[i] = zero_
                x_rho_w[i] = (metric_factor[column, i + 1, layer] *
                              derivative_hat[i + 1, node] *
                              dp_drho_theta[column, node, layer])
            end
        else
            # Column of `rho_w`: `b_rho = A13[:, node]`, whose diagonal entry vanishes at
            # the interfaces between elements since `S / 2` has a zero diagonal,
            # `b_rho_theta = A23[:, node]`, and `b_rho_w = 0`
            for i in 1:NNODES
                x_rho[i] = ((is_interface && i == node) ? zero_ :
                            metric_factor[column, i, layer] * derivative_hat[i, node])
            end
            for i in 1:(NNODES - 2)
                x_rho_theta[i] = (metric_factor[column, i + 1, layer] *
                                  derivative_hat[i + 1, node] *
                                  theta[column, node, layer])
            end
            for i in 1:(NNODES - 2)
                x_rho_w[i] = gravity[column, i + 1, layer] * x_rho[i + 1]
                for j in 1:(NNODES - 2)
                    x_rho_w[i] += (metric_factor[column, i + 1, layer] *
                                   derivative_hat[i + 1, j + 1] *
                                   dp_drho_theta[column, j + 1, layer]) * x_rho_theta[j]
                end
                x_rho_w[i] *= -gamma_dt
            end
        end

        schur_complement_boundary_unknown!(schur_band, x_rho, x_rho_theta, x_rho_w,
                                           schur_element, gravity, theta, dp_drho_theta,
                                           metric_factor, derivative_hat, gamma_dt,
                                           Val(NNODES), index_local, offset_boundary,
                                           column, layer)
    end

    # Gravity couples the density in the edge nodes to the boundary unknowns of the
    # neighboring elements through the density entries of `B` at the lower interface (the
    # layer below as `ll` and this layer as `rr`) and the upper interface (this layer as
    # `ll` and the layer above as `rr`)
    gravity_1 = gravity[column, 1, layer]
    gravity_M = gravity[column, NNODES, layer]
    if layer > 1
        dp_drho_theta_ll = dp_drho_theta[column, NNODES, layer - 1]
        factor_rr = 0.5f0 * metric_factor[column, 1, layer] * inv_boundary_weight
        schur_band[column, 7, offset_boundary - 1] += gravity_1 * gamma_dt * factor_rr
        schur_band[column, 6, offset_boundary] += gravity_1 * gamma_dt *
                                                  dp_drho_theta_ll *
                                                  inv_speed_of_sound *
                                                  factor_rr
    end
    if layer < nlayers_
        dp_drho_theta_rr = dp_drho_theta[column, 1, layer + 1]
        factor_ll = 0.5f0 * metric_factor[column, NNODES, layer] * inv_boundary_weight
        schur_band[column, 2, offset_boundary + 5] += gravity_M * gamma_dt *
                                                      dp_drho_theta_rr *
                                                      inv_speed_of_sound * factor_ll
        schur_band[column, 1, offset_boundary + 6] += -gravity_M * gamma_dt * factor_ll
    end
    return nothing
end

# Add `-C A^{-1} b` of the column `b` of `B` of the boundary unknown `index_local` of
# element `layer` to `schur_band`. With the blocks `A_ij` of the element block `A_e` of the
# reference in the docstring of `VerticalJacobianTwoSided`, `x = A^{-1} b` follows from the
# local Schur complement `S_e`:
#   S_e x_rho_w = b_rho_w - A31 A11^{-1} b_rho - A32 A22^{-1} b_rho_theta,
#   x_rho = A11^{-1} (b_rho - A13 x_rho_w),
#   x_rho_theta = A22^{-1} (b_rho_theta - A23 x_rho_w),
# with `A11 = A22 = I / gamma_dt`. On entry, `x_rho` and `x_rho_theta` hold the components
# of `b` and `x_rho_w` the right-hand side of `S_e x_rho_w`, assembled by the caller
# (`schur_complement_boundary_per_element!`), which skips the terms that vanish for the columns
# of `rho_theta`. `x_rho_theta` and `x_rho_w` are overwritten by the components of `x`. Only
# the density in the edge nodes is needed, since `C` couples the density only there.
Base.@propagate_inbounds function schur_complement_boundary_unknown!(schur_band, x_rho,
                                                                     x_rho_theta,
                                                                     x_rho_w,
                                                                     schur_element,
                                                                     gravity, theta,
                                                                     dp_drho_theta,
                                                                     metric_factor,
                                                                     derivative_hat,
                                                                     gamma_dt,
                                                                     ::Val{NNODES},
                                                                     index_local,
                                                                     offset_boundary,
                                                                     column,
                                                                     layer) where {NNODES}
    RealT = eltype(schur_band)
    zero_ = zero(RealT)
    metric_factor_1 = metric_factor[column, 1, layer]
    metric_factor_M = metric_factor[column, NNODES, layer]

    # Solve `S_e x_rho_w` with the LU factors of `S_e` in place
    ldiv_dense!(x_rho_w, schur_element, Val(NNODES - 2), column, layer)

    # Density in the edge nodes
    x_rho_1 = x_rho[1]
    x_rho_M = x_rho[NNODES]
    for j in 1:(NNODES - 2)
        x_rho_1 += -(metric_factor_1 * derivative_hat[1, j + 1]) * x_rho_w[j]
        x_rho_M += -(metric_factor_M * derivative_hat[NNODES, j + 1]) * x_rho_w[j]
    end
    x_rho_1 *= gamma_dt
    x_rho_M *= gamma_dt

    # Interior potential temperature density
    for i in 1:(NNODES - 2)
        for j in 1:(NNODES - 2)
            x_rho_theta[i] += -(metric_factor[column, i + 1, layer] *
                                derivative_hat[i + 1, j + 1] *
                                theta[column, j + 1, layer]) * x_rho_w[j]
        end
        x_rho_theta[i] *= gamma_dt
    end

    # Entries `-C x` in the rows `rho_theta_1`, `rho_w_1`, `rho_w_M`, and `rho_theta_M`
    schur_rho_theta_1 = zero_
    schur_rho_w_1 = zero_
    schur_rho_w_M = zero_
    schur_rho_theta_M = zero_
    for j in 1:(NNODES - 2)
        theta_j = theta[column, j + 1, layer]
        dp_drho_theta_j = dp_drho_theta[column, j + 1, layer]
        schur_rho_theta_1 += -(metric_factor_1 * derivative_hat[1, j + 1] * theta_j) *
                             x_rho_w[j]
        schur_rho_theta_M += -(metric_factor_M * derivative_hat[NNODES, j + 1] *
                               theta_j) * x_rho_w[j]
        schur_rho_w_1 += -(metric_factor_1 * derivative_hat[1, j + 1] *
                           dp_drho_theta_j) * x_rho_theta[j]
        schur_rho_w_M += -(metric_factor_M * derivative_hat[NNODES, j + 1] *
                           dp_drho_theta_j) * x_rho_theta[j]
    end
    schur_rho_w_1 -= gravity[column, 1, layer] * x_rho_1
    schur_rho_w_M -= gravity[column, NNODES, layer] * x_rho_M

    # The entry `(offset_boundary + i, offset_boundary + index_local)` is stored at
    # `schur_band[column, i - index_local + 4, offset_boundary + index_local]`
    index_column = offset_boundary + index_local
    schur_band[column, 5 - index_local, index_column] += schur_rho_theta_1
    schur_band[column, 6 - index_local, index_column] += schur_rho_w_1
    schur_band[column, 7 - index_local, index_column] += schur_rho_w_M
    schur_band[column, 8 - index_local, index_column] += schur_rho_theta_M
    return nothing
end

@kernel inbounds=true function schur_complement_boundary_KAkernel!(schur_band,
                                                                   @Const(schur_element),
                                                                   @Const(gravity),
                                                                   @Const(theta),
                                                                   @Const(dp_drho_theta),
                                                                   @Const(metric_factor),
                                                                   @Const(derivative_hat),
                                                                   gamma_dt,
                                                                   inv_speed_of_sound,
                                                                   inv_boundary_weight,
                                                                   nlayers_,
                                                                   ::Val{NNODES}) where {NNODES}
    column = @index(Global)
    RealT = eltype(schur_band)
    x_rho = MVector{NNODES, RealT}(undef)
    x_rho_theta = MVector{NNODES - 2, RealT}(undef)
    x_rho_w = MVector{NNODES - 2, RealT}(undef)

    for layer in 1:nlayers_
        schur_complement_boundary_per_element!(schur_band, x_rho, x_rho_theta, x_rho_w,
                                               schur_element, gravity, theta,
                                               dp_drho_theta, metric_factor,
                                               derivative_hat, gamma_dt,
                                               inv_speed_of_sound, inv_boundary_weight,
                                               nlayers_, Val(NNODES), column, layer)
    end
end

function schur_complement_boundary!(backend::Backend,
                                    jacobian::VerticalJacobianTwoSidedCache,
                                    gamma_dt)
    (; schur_element, ncolumns, schur_band, gravity, theta, dp_drho_theta, metric_factor,
    derivative_hat, inv_boundary_weight) = jacobian
    nnodes_ = nnodes(jacobian)
    nlayers_ = nlayers(jacobian)
    inv_speed_of_sound = inv(jacobian.speed_of_sound)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck check_axes_columns(jacobian)

    kernel! = schur_complement_boundary_KAkernel!(backend)
    kernel!(schur_band, schur_element, gravity, theta, dp_drho_theta, metric_factor,
            derivative_hat, gamma_dt, inv_speed_of_sound, inv_boundary_weight, nlayers_,
            Val(nnodes_);
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Banded LU factorization of the boundary Schur complement and the corresponding solve.
# `schur_band[column, :, :]` stores the 3 lower and 3 upper diagonals in the usual banded
# layout, i.e., the entry `(i, j)` at `schur_band[column, i - j + 4, j]`.

# The diagonal is stored in the row `bandwidth_upper + 1` of the banded layout, and the
# subdiagonals below it.
Base.@propagate_inbounds function lu_band_per_column!(A, bandwidth_lower,
                                                      bandwidth_upper,
                                                      ::Val{N}, column) where {N}
    diagonal = bandwidth_upper + 1
    first_subdiagonal = diagonal + 1
    last_subdiagonal = diagonal + bandwidth_lower
    for k in 1:(N - 1)
        # One division and `bandwidth_lower` multiplications instead of `bandwidth_lower`
        # divisions
        inv_kk = inv(A[column, diagonal, k])
        for i in first_subdiagonal:last_subdiagonal
            A[column, i, k] *= inv_kk
        end
        for i in first_subdiagonal:last_subdiagonal
            for j in (k + 1):min(k + bandwidth_upper, N)
                A[column, i + k - j, j] -= A[column, i, k] *
                                           A[column, k + diagonal - j, j]
            end
        end
    end
    return nothing
end

@kernel inbounds=true function lu_band_KAkernel!(schur_band, ::Val{N}) where {N}
    column = @index(Global)

    lu_band_per_column!(schur_band, 3, 3, Val(N), column)
end

function lu_band!(backend::Backend, jacobian::VerticalJacobianTwoSidedCache)
    (; ncolumns, schur_band) = jacobian
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck check_axes_columns(jacobian)
    kernel! = lu_band_KAkernel!(backend)
    kernel!(schur_band, Val(4 * nlayers(jacobian));
            ndrange = ncolumns, column_kernel_workgroupsize()...)
    return nothing
end

# Solve with the banded LU factors `A = schur_band` and the right-hand side `b = schur_rhs`
Base.@propagate_inbounds function ldiv_band_per_column!(b, A, bandwidth_lower,
                                                        bandwidth_upper, ::Val{N},
                                                        column) where {N}
    diagonal = bandwidth_upper + 1
    # Ly = b
    for k in 1:(N - 1)
        for i in (k + 1):min(k + bandwidth_lower, N)
            b[column, i] -= A[column, i - k + diagonal, k] * b[column, k]
        end
    end
    # Ux = y
    for k in N:-1:1
        b[column, k] /= A[column, diagonal, k]
        for i in max(k - bandwidth_upper, 1):(k - 1)
            b[column, i] -= A[column, diagonal + i - k, k] * b[column, k]
        end
    end
    return nothing
end

@kernel inbounds=true function ldiv_band_KAkernel!(schur_rhs, @Const(schur_band),
                                                   ::Val{N}) where {N}
    column = @index(Global)

    ldiv_band_per_column!(schur_rhs, schur_band, 3, 3, Val(N), column)
end

function ldiv_band!(backend::Backend, jacobian::VerticalJacobianTwoSidedCache)
    (; ncolumns, schur_band, schur_rhs) = jacobian
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck check_axes_columns(jacobian)
    kernel! = ldiv_band_KAkernel!(backend)
    kernel!(schur_rhs, schur_band, Val(4 * nlayers(jacobian));
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
Base.@propagate_inbounds function ldiv_interior_forward_per_element!(schur_rhs, x_rho,
                                                                     x_rho_theta,
                                                                     x_rho_w,
                                                                     schur_rhs_element,
                                                                     b, schur_element,
                                                                     gravity, theta,
                                                                     dp_drho_theta,
                                                                     metric_factor,
                                                                     derivative_hat,
                                                                     gamma_dt,
                                                                     ::Val{NNODES},
                                                                     column,
                                                                     layer) where {NNODES}
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
                           dp_drho_theta[column, j + 1, layer]) * x_rho_theta[j]
        end
    end
    ldiv_dense!(x_rho_w, schur_element, Val(NNODES - 2), column, layer)
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
        dp_drho_theta_j = dp_drho_theta[column, j + 1, layer]
        schur_rhs_element[1] += -(metric_factor[column, 1, layer] *
                                  derivative_hat[1, j + 1] * theta_j) *
                                x_rho_w[j]
        schur_rhs_element[4] += -(metric_factor[column, NNODES, layer] *
                                  derivative_hat[NNODES, j + 1] *
                                  theta_j) * x_rho_w[j]
        schur_rhs_element[2] += -(metric_factor[column, 1, layer] *
                                  derivative_hat[1, j + 1] * dp_drho_theta_j) *
                                x_rho_theta[j]
        schur_rhs_element[3] += -(metric_factor[column, NNODES, layer] *
                                  derivative_hat[NNODES, j + 1] *
                                  dp_drho_theta_j) * x_rho_theta[j]
    end
    schur_rhs_element[2] += -gravity[column, 1, layer] * x_rho_1
    schur_rhs_element[3] += -gravity[column, NNODES, layer] * x_rho_M
    schur_rhs[column, offset_boundary + 1] = schur_rhs_element[1]
    schur_rhs[column, offset_boundary + 2] = schur_rhs_element[2]
    schur_rhs[column, offset_boundary + 3] = schur_rhs_element[3]
    schur_rhs[column, offset_boundary + 4] = schur_rhs_element[4]
    return nothing
end

@kernel inbounds=true function ldiv_interior_forward_KAkernel!(schur_rhs, @Const(b),
                                                               @Const(schur_element),
                                                               @Const(gravity),
                                                               @Const(theta),
                                                               @Const(dp_drho_theta),
                                                               @Const(metric_factor),
                                                               @Const(derivative_hat),
                                                               gamma_dt,
                                                               ::Val{NNODES}) where {NNODES}
    column, layer = @index(Global, NTuple)
    x_rho = MVector{NNODES, eltype(b)}(undef)
    x_rho_theta = MVector{NNODES - 2, eltype(b)}(undef)
    x_rho_w = MVector{NNODES - 2, eltype(b)}(undef)
    schur_rhs_element = MVector{4, eltype(b)}(undef)

    ldiv_interior_forward_per_element!(schur_rhs, x_rho, x_rho_theta, x_rho_w,
                                       schur_rhs_element, b, schur_element, gravity,
                                       theta, dp_drho_theta, metric_factor,
                                       derivative_hat, gamma_dt, Val(NNODES), column,
                                       layer)
end

function ldiv_interior_forward!(backend::Backend,
                                jacobian::VerticalJacobianTwoSidedCache, b,
                                gamma_dt)
    (; schur_rhs, schur_element, ncolumns, theta, dp_drho_theta, metric_factor, derivative_hat,
    gravity) = jacobian
    nnodes_ = nnodes(jacobian)
    nlayers_ = nlayers(jacobian)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        check_axes_columns(jacobian)
        Trixi.check_axes(b, axes(jacobian.u_vertical))
    end

    kernel! = ldiv_interior_forward_KAkernel!(backend)
    kernel!(schur_rhs, b, schur_element, gravity, theta, dp_drho_theta, metric_factor,
            derivative_hat, gamma_dt, Val(nnodes_);
            ndrange = (ncolumns, nlayers_), layer_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Backward substitution: the solution of the boundary degrees of freedom `schur_rhs` and
# of the interior degrees of freedom

# Solution of element `layer` of the column `column`; `x_rho`, `x_rho_theta`, and `x_rho_w`
# hold the interior components of the three variables, and `rho_theta_1`, `rho_w_1`,
# `rho_w_M`, and `rho_theta_M` the boundary components from `schur_rhs`
Base.@propagate_inbounds function ldiv_interior_backward_per_element!(b, x_rho,
                                                                      x_rho_theta,
                                                                      x_rho_w,
                                                                      schur_rhs,
                                                                      schur_element,
                                                                      gravity, theta,
                                                                      dp_drho_theta,
                                                                      metric_factor,
                                                                      derivative_hat,
                                                                      gamma_dt,
                                                                      inv_speed_of_sound,
                                                                      inv_boundary_weight,
                                                                      nlayers_,
                                                                      ::Val{NNODES},
                                                                      column,
                                                                      layer) where {NNODES}
    RealT = eltype(b)
    zero_ = zero(RealT)
    # Offset of the 4 boundary unknowns of this layer in the condensed system
    offset_boundary = (layer - 1) * 4
    theta_1 = theta[column, 1, layer]
    theta_M = theta[column, NNODES, layer]
    dp_drho_theta_1 = dp_drho_theta[column, 1, layer]
    dp_drho_theta_M = dp_drho_theta[column, NNODES, layer]
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
    # Density entries of `B` from the LMARS mass flux at the lower interface of the element,
    # with the layer below as `ll` and this layer as `rr`; it vanishes at the wall
    if layer > 1
        factor_rr = 0.5f0 * metric_factor[column, 1, layer] * inv_boundary_weight
        x_rho[1] -= dp_drho_theta_1 * inv_speed_of_sound * factor_rr *
                    schur_rhs[column, offset_boundary + 1]
    end
    # Density entries of `B` from the LMARS mass flux at the upper interface of the element,
    # with this layer as `ll` and the layer above as `rr`; it vanishes at the wall
    if layer < nlayers_
        factor_ll = 0.5f0 * metric_factor[column, NNODES, layer] * inv_boundary_weight
        x_rho[NNODES] -= dp_drho_theta_M * inv_speed_of_sound * factor_ll *
                         schur_rhs[column, offset_boundary + 4]
    end
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
        if layer < nlayers_ && i == NNODES
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
                       derivative_hat[i + 1, 1] * dp_drho_theta_1) *
                      rho_theta_1 +
                      (metric_factor[column, i + 1, layer] *
                       derivative_hat[i + 1, NNODES] *
                       dp_drho_theta_M) * rho_theta_M
    end
    # Density entries of `B` coupling to the boundary unknowns of the neighboring elements
    # at the lower interface (the layer below as `ll` and this layer as `rr`) and the upper
    # interface (this layer as `ll` and the layer above as `rr`)
    if layer > 1
        dp_drho_theta_ll = dp_drho_theta[column, NNODES, layer - 1]
        factor_rr = 0.5f0 * metric_factor[column, 1, layer] * inv_boundary_weight
        x_rho[1] += factor_rr * schur_rhs[column, offset_boundary - 1] +
                    dp_drho_theta_ll * inv_speed_of_sound * factor_rr *
                    schur_rhs[column, offset_boundary]
    end
    if layer < nlayers_
        dp_drho_theta_rr = dp_drho_theta[column, 1, layer + 1]
        factor_ll = 0.5f0 * metric_factor[column, NNODES, layer] * inv_boundary_weight
        x_rho[NNODES] += dp_drho_theta_rr * inv_speed_of_sound * factor_ll *
                         schur_rhs[column, offset_boundary + 5] -
                         factor_ll * schur_rhs[column, offset_boundary + 6]
    end
    for i in 1:(NNODES - 2)
        x_rho_w[i] += -gamma_dt * gravity[column, i + 1, layer] * x_rho[i + 1]
        for j in 1:(NNODES - 2)
            x_rho_w[i] += -gamma_dt *
                          (metric_factor[column, i + 1, layer] *
                           derivative_hat[i + 1, j + 1] *
                           dp_drho_theta[column, j + 1, layer]) * x_rho_theta[j]
        end
    end
    ldiv_dense!(x_rho_w, schur_element, Val(NNODES - 2), column, layer)
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

@kernel inbounds=true function ldiv_interior_backward_KAkernel!(b, @Const(schur_rhs),
                                                                @Const(schur_element),
                                                                @Const(gravity),
                                                                @Const(theta),
                                                                @Const(dp_drho_theta),
                                                                @Const(metric_factor),
                                                                @Const(derivative_hat),
                                                                gamma_dt,
                                                                inv_speed_of_sound,
                                                                inv_boundary_weight,
                                                                nlayers_,
                                                                ::Val{NNODES}) where {NNODES}
    column, layer = @index(Global, NTuple)
    x_rho = MVector{NNODES, eltype(b)}(undef)
    x_rho_theta = MVector{NNODES - 2, eltype(b)}(undef)
    x_rho_w = MVector{NNODES - 2, eltype(b)}(undef)

    ldiv_interior_backward_per_element!(b, x_rho, x_rho_theta, x_rho_w, schur_rhs,
                                        schur_element, gravity, theta, dp_drho_theta,
                                        metric_factor, derivative_hat, gamma_dt,
                                        inv_speed_of_sound, inv_boundary_weight,
                                        nlayers_, Val(NNODES), column, layer)
end

function ldiv_interior_backward!(backend::Backend,
                                 jacobian::VerticalJacobianTwoSidedCache,
                                 b,
                                 gamma_dt)
    (; schur_element, schur_rhs, ncolumns, gravity, theta, dp_drho_theta, metric_factor,
    derivative_hat, inv_boundary_weight) = jacobian
    nnodes_ = nnodes(jacobian)
    nlayers_ = nlayers(jacobian)
    inv_speed_of_sound = inv(jacobian.speed_of_sound)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        check_axes_columns(jacobian)
        Trixi.check_axes(b, axes(jacobian.u_vertical))
    end

    kernel! = ldiv_interior_backward_KAkernel!(backend)
    kernel!(b, schur_rhs, schur_element, gravity, theta, dp_drho_theta, metric_factor,
            derivative_hat, gamma_dt, inv_speed_of_sound, inv_boundary_weight, nlayers_,
            Val(nnodes_);
            ndrange = (ncolumns, nlayers_), layer_kernel_workgroupsize()...)
    return nothing
end

###############################################################################
# Dense LU factorization and solve of the local Schur complement `S_e` of one element
# (indexed `[column, i, j, layer]`, `N = nnodes - 2`)
Base.@propagate_inbounds function lu_dense!(A, ::Val{N}, column, layer) where {N}
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

Base.@propagate_inbounds function ldiv_dense!(b, A, ::Val{N}, column, layer) where {N}
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
