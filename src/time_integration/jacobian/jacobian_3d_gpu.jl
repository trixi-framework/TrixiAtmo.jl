# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Copies of the data of one node between Trixi.jl's layout and the columns for the
# `CompressibleEulerPotentialTemperatureEquationsWithGravity3D`, called by the copy kernels
# in `jacobian_gpu.jl`

# Copy the vertical variables of the node `(i, j, k, element)` of `u`, multiplied by
# `scaling`, where the vertical momentum is the component of the momentum in the
# `vertical_direction` `nu` (see `compute_vertical_frame`). The `scaling` is one for the
# state and the normalized norm `n` for the right-hand side of the linear system, i.e.,
# the diagonal similarity transformation `N` of the vertical operator.
Base.@propagate_inbounds function copy_node_to_columns!(u_vertical, u, scaling,
                                                        vertical_frame,
                                                        horizontal_index_of_element,
                                                        vertical_index_of_element,
                                                        nnodes_,
                                                        i, j, k, element)
    column = (horizontal_index_of_element[element] - 1) * nnodes_ * nnodes_ +
             (j - 1) * nnodes_ + i
    layer = vertical_index_of_element[element]
    vertical_direction = SVector(vertical_frame[1, i, j, k, element],
                                 vertical_frame[2, i, j, k, element],
                                 vertical_frame[3, i, j, k, element])
    rho_v = SVector(u[2, i, j, k, element], u[3, i, j, k, element],
                    u[4, i, j, k, element])
    u_vertical[column, 1, k, layer] = scaling * u[1, i, j, k, element]
    u_vertical[column, 2, k, layer] = scaling * dot(vertical_direction, rho_v)
    u_vertical[column, 3, k, layer] = scaling * u[5, i, j, k, element]
    return nothing
end

# Copy the solution of the linear system, divided by the normalized norm `n`, to `b`. The
# horizontal components of the momentum are multiplied by `gamma_dt`, i.e., solved with
# `J = 0`, and the component in the `vertical_direction` is replaced by the solution.
Base.@propagate_inbounds function copy_node_from_columns!(b, u_vertical, vertical_frame,
                                                          horizontal_index_of_element,
                                                          vertical_index_of_element,
                                                          nnodes_, gamma_dt,
                                                          i, j, k, element)
    column = (horizontal_index_of_element[element] - 1) * nnodes_ * nnodes_ +
             (j - 1) * nnodes_ + i
    layer = vertical_index_of_element[element]
    vertical_direction = SVector(vertical_frame[1, i, j, k, element],
                                 vertical_frame[2, i, j, k, element],
                                 vertical_frame[3, i, j, k, element])
    inv_area = inv(vertical_frame[4, i, j, k, element])
    rho_v = SVector(b[2, i, j, k, element], b[3, i, j, k, element],
                    b[4, i, j, k, element])
    rho_w = inv_area * u_vertical[column, 2, k, layer]
    rho_v = gamma_dt * rho_v +
            (rho_w - gamma_dt * dot(vertical_direction, rho_v)) * vertical_direction
    b[1, i, j, k, element] = inv_area * u_vertical[column, 1, k, layer]
    b[2, i, j, k, element] = rho_v[1]
    b[3, i, j, k, element] = rho_v[2]
    b[4, i, j, k, element] = rho_v[3]
    b[5, i, j, k, element] = inv_area * u_vertical[column, 3, k, layer]
    return nothing
end
end # @muladd
