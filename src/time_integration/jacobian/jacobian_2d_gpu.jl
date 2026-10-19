# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Copies of the data of one node between Trixi.jl's layout and the columns for the
# `CompressibleEulerPotentialTemperatureEquationsWithGravity2D`, called by the copy kernels
# in `jacobian_gpu.jl` (see `jacobian_3d_gpu.jl`)

# In 2D, the vertical node index is `j`, and the variables are
# `(rho, rho_v1, rho_v2, rho_theta, phi)`
Base.@propagate_inbounds function copy_node_to_columns!(u_vertical, u, scaling,
                                                        vertical_frame,
                                                        horizontal_index_of_element,
                                                        vertical_index_of_element,
                                                        nnodes_,
                                                        i, j, element)
    column = (horizontal_index_of_element[element] - 1) * nnodes_ + i
    layer = vertical_index_of_element[element]
    vertical_direction = SVector(vertical_frame[1, i, j, element],
                                 vertical_frame[2, i, j, element])
    rho_v = SVector(u[2, i, j, element], u[3, i, j, element])
    u_vertical[column, 1, j, layer] = scaling * u[1, i, j, element]
    u_vertical[column, 2, j, layer] = scaling * dot(vertical_direction, rho_v)
    u_vertical[column, 3, j, layer] = scaling * u[4, i, j, element]
    return nothing
end

Base.@propagate_inbounds function copy_node_from_columns!(b, u_vertical, vertical_frame,
                                                          horizontal_index_of_element,
                                                          vertical_index_of_element,
                                                          nnodes_, gamma_dt,
                                                          i, j, element)
    column = (horizontal_index_of_element[element] - 1) * nnodes_ + i
    layer = vertical_index_of_element[element]
    vertical_direction = SVector(vertical_frame[1, i, j, element],
                                 vertical_frame[2, i, j, element])
    inv_area = inv(vertical_frame[3, i, j, element])
    rho_v = SVector(b[2, i, j, element], b[3, i, j, element])
    rho_w = inv_area * u_vertical[column, 2, j, layer]
    rho_v = gamma_dt * rho_v +
            (rho_w - gamma_dt * dot(vertical_direction, rho_v)) * vertical_direction
    b[1, i, j, element] = inv_area * u_vertical[column, 1, j, layer]
    b[2, i, j, element] = rho_v[1]
    b[3, i, j, element] = rho_v[2]
    b[4, i, j, element] = inv_area * u_vertical[column, 3, j, layer]
    return nothing
end
end # @muladd
