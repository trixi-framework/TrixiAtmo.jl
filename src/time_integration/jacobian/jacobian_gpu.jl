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
# momentum), and `3 => rho_theta`. The copies of the data of one node are in
# `jacobian_3d_gpu.jl`. The kernels launch one thread per node, i.e.,
# `indices = (i, j, k, element)`.

@kernel inbounds=true function copy_state_to_columns_KAkernel!(u_vertical, @Const(u),
                                                               @Const(vertical_frame),
                                                               @Const(horizontal_index_of_element),
                                                               @Const(vertical_index_of_element),
                                                               nnodes_)
    indices = @index(Global, NTuple)
    copy_node_to_columns!(u_vertical, u, one(eltype(u_vertical)), vertical_frame,
                          horizontal_index_of_element, vertical_index_of_element,
                          nnodes_, indices...)
end

@kernel inbounds=true function copy_rhs_to_columns_KAkernel!(u_vertical, @Const(b),
                                                             @Const(vertical_frame),
                                                             @Const(horizontal_index_of_element),
                                                             @Const(vertical_index_of_element),
                                                             nnodes_)
    indices = @index(Global, NTuple)
    # The last slot of the vertical frame holds the normalized norm `n`
    copy_node_to_columns!(u_vertical, b, vertical_frame[end, indices...],
                          vertical_frame, horizontal_index_of_element,
                          vertical_index_of_element, nnodes_, indices...)
end

@kernel inbounds=true function copy_solution_from_columns_KAkernel!(b,
                                                                    @Const(u_vertical),
                                                                    @Const(vertical_frame),
                                                                    @Const(horizontal_index_of_element),
                                                                    @Const(vertical_index_of_element),
                                                                    nnodes_, gamma_dt)
    indices = @index(Global, NTuple)
    copy_node_from_columns!(b, u_vertical, vertical_frame, horizontal_index_of_element,
                            vertical_index_of_element, nnodes_, gamma_dt, indices...)
end

# Copy the state `u`, from which the Jacobian is assembled, to the columns
function copy_state_to_columns!(backend::Backend, u_vertical, u, jacobian, semi)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck Trixi.check_axes(u, mesh, equations, solver, cache)
    (; vertical_frame, horizontal_index_of_element, vertical_index_of_element) = jacobian
    nnodes_ = number_of_nodes(jacobian)
    kernel! = copy_state_to_columns_KAkernel!(backend)
    kernel!(u_vertical, u, vertical_frame, horizontal_index_of_element,
            vertical_index_of_element, nnodes_;
            ndrange = Base.tail(size(u)))
    return nothing
end

# Copy the right-hand side `b` of the linear system, multiplied by the normalized norm
# `n`, to the columns
function copy_rhs_to_columns!(backend::Backend, u_vertical, b, jacobian, semi)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck Trixi.check_axes(b, mesh, equations, solver, cache)
    (; vertical_frame, horizontal_index_of_element, vertical_index_of_element) = jacobian
    nnodes_ = number_of_nodes(jacobian)
    kernel! = copy_rhs_to_columns_KAkernel!(backend)
    kernel!(u_vertical, b, vertical_frame, horizontal_index_of_element,
            vertical_index_of_element, nnodes_;
            ndrange = Base.tail(size(b)))
    return nothing
end

# Copy the solution of the linear system, divided by the normalized norm `n`, from the
# columns to `b`
function copy_solution_from_columns!(backend::Backend, b, u_vertical, jacobian, semi,
                                     gamma_dt)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck Trixi.check_axes(b, mesh, equations, solver, cache)
    (; vertical_frame, horizontal_index_of_element, vertical_index_of_element) = jacobian
    nnodes_ = number_of_nodes(jacobian)
    kernel! = copy_solution_from_columns_KAkernel!(backend)
    kernel!(b, u_vertical, vertical_frame, horizontal_index_of_element,
            vertical_index_of_element, nnodes_, gamma_dt;
            ndrange = Base.tail(size(b)))
    return nothing
end
end # @muladd
