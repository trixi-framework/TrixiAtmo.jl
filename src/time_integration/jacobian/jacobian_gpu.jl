# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Launch parameters and copy kernels shared by all vertical Jacobians

# Launch keywords of the column-parallel kernels (one thread per vertical column) and of the
# layer-parallel kernels (one thread per `(column, layer)`): the workgroup sizes measured
# fastest on an NVIDIA GH200 and on an AMD MI210.
@inline column_kernel_workgroupsize() = (workgroupsize = 64,)
@inline layer_kernel_workgroupsize() = (workgroupsize = 128,)

# Copy data between Trixi.jl's layout `[v, i, j, k, element]` and the layout of the columns
# `[column, v, node, layer]`, with the variables `1 => rho`, `2 => rho_w` (vertical
# momentum), and `3 => rho_theta`, via `cons2hevi` and `hevi2cons` of the equations (see
# `jacobian_3d_gpu.jl`). The kernels launch one thread per node `(i, j, k, element)`.

# Check whether the arrays in `jacobian` have the axes we assume they must have in the copy
# kernels between Trixi.jl's layout and the columns.
function check_axes_element_layout(jacobian, solver, cache)
    (; vertical_direction_and_norm, horizontal_index_of_element,
    vertical_index_of_element) = jacobian
    nodes = eachnode(solver)
    elements = eachelement(solver, cache)
    Trixi.check_axes(vertical_direction_and_norm,
                     (Base.OneTo(4), nodes, nodes, nodes, elements))
    Trixi.check_axes(horizontal_index_of_element, (elements,))
    Trixi.check_axes(vertical_index_of_element, (elements,))
    check_axes_columns(jacobian)
    return nothing
end

# Copy the state `u`, from which the Jacobian is assembled, to the columns
function copy_state_to_columns!(backend::Backend, u_vertical, u, jacobian, semi)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        Trixi.check_axes(u, mesh, equations, solver, cache)
        check_axes_element_layout(jacobian, solver, cache)
        Trixi.check_axes(u_vertical, axes(jacobian.u_vertical))
    end
    (; vertical_direction_and_norm, horizontal_index_of_element,
    vertical_index_of_element) = jacobian
    nnodes_ = nnodes(jacobian)
    kernel! = copy_state_to_columns_KAkernel!(backend)
    kernel!(u_vertical, u, vertical_direction_and_norm, horizontal_index_of_element,
            vertical_index_of_element, nnodes_, equations, solver;
            ndrange = Base.tail(size(u)))
    return nothing
end

@kernel inbounds=true function copy_state_to_columns_KAkernel!(u_vertical, @Const(u),
                                                               @Const(vertical_direction_and_norm),
                                                               @Const(horizontal_index_of_element),
                                                               @Const(vertical_index_of_element),
                                                               nnodes_, equations,
                                                               solver)
    i, j, k, element = @index(Global, NTuple)
    column, layer = get_column_and_layer(horizontal_index_of_element,
                                         vertical_index_of_element, nnodes_, i, j,
                                         element)
    vertical_direction = get_vertical_direction(vertical_direction_and_norm,
                                                i, j, k, element)
    u_node = get_node_vars(u, equations, solver, i, j, k, element)
    set_column_node_vars!(u_vertical, cons2hevi(u_node, vertical_direction, equations),
                          column, k, layer)
end

# Copy the right-hand side `b` of the linear system, multiplied by the norm
# `n`, to the columns
function copy_rhs_to_columns!(backend::Backend, u_vertical, b, jacobian, semi)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        Trixi.check_axes(b, mesh, equations, solver, cache)
        check_axes_element_layout(jacobian, solver, cache)
        Trixi.check_axes(u_vertical, axes(jacobian.u_vertical))
    end
    (; vertical_direction_and_norm, horizontal_index_of_element,
    vertical_index_of_element) = jacobian
    nnodes_ = nnodes(jacobian)
    kernel! = copy_rhs_to_columns_KAkernel!(backend)
    kernel!(u_vertical, b, vertical_direction_and_norm, horizontal_index_of_element,
            vertical_index_of_element, nnodes_, equations, solver;
            ndrange = Base.tail(size(b)))
    return nothing
end

@kernel inbounds=true function copy_rhs_to_columns_KAkernel!(u_vertical, @Const(b),
                                                             @Const(vertical_direction_and_norm),
                                                             @Const(horizontal_index_of_element),
                                                             @Const(vertical_index_of_element),
                                                             nnodes_, equations, solver)
    i, j, k, element = @index(Global, NTuple)
    column, layer = get_column_and_layer(horizontal_index_of_element,
                                         vertical_index_of_element, nnodes_, i, j,
                                         element)
    vertical_direction = get_vertical_direction(vertical_direction_and_norm,
                                                i, j, k, element)
    norm_Ja3 = get_vertical_norm(vertical_direction_and_norm, i, j, k, element)
    b_node = get_node_vars(b, equations, solver, i, j, k, element)
    set_column_node_vars!(u_vertical,
                          norm_Ja3 * cons2hevi(b_node, vertical_direction, equations),
                          column, k, layer)
end

# Copy the solution of the linear system, divided by the norm `n`, from the
# columns to `b`
function copy_solution_from_columns!(backend::Backend, b, u_vertical, jacobian, semi,
                                     gamma_dt)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        Trixi.check_axes(b, mesh, equations, solver, cache)
        check_axes_element_layout(jacobian, solver, cache)
        Trixi.check_axes(u_vertical, axes(jacobian.u_vertical))
    end
    (; vertical_direction_and_norm, horizontal_index_of_element,
    vertical_index_of_element) = jacobian
    nnodes_ = nnodes(jacobian)
    kernel! = copy_solution_from_columns_KAkernel!(backend)
    kernel!(b, u_vertical, vertical_direction_and_norm, horizontal_index_of_element,
            vertical_index_of_element, nnodes_, gamma_dt, equations, solver;
            ndrange = Base.tail(size(b)))
    return nothing
end

@kernel inbounds=true function copy_solution_from_columns_KAkernel!(b,
                                                                    @Const(u_vertical),
                                                                    @Const(vertical_direction_and_norm),
                                                                    @Const(horizontal_index_of_element),
                                                                    @Const(vertical_index_of_element),
                                                                    nnodes_, gamma_dt,
                                                                    equations, solver)
    i, j, k, element = @index(Global, NTuple)
    column, layer = get_column_and_layer(horizontal_index_of_element,
                                         vertical_index_of_element, nnodes_, i, j,
                                         element)
    vertical_direction = get_vertical_direction(vertical_direction_and_norm,
                                                i, j, k, element)
    inv_norm_Ja3 = inv(get_vertical_norm(vertical_direction_and_norm, i, j, k, element))
    x_vertical = inv_norm_Ja3 * get_column_node_vars(u_vertical, column, k, layer)
    b_node = get_node_vars(b, equations, solver, i, j, k, element)
    set_node_vars!(b,
                   hevi2cons(x_vertical, b_node, vertical_direction, gamma_dt,
                             equations),
                   equations, solver, i, j, k, element)
end
end # @muladd
