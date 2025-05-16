using Trixi: get_node_vars, @batch, set_node_vars!, get_node_coords, each_quad_node
using LinearAlgebra
using StaticArrays

@muladd begin
    function rain_limiter_dg2d!(u, equations::AbstractCompressibleRainyEulerEquations,
                                dg::DGSEM, cache, mesh)

        # iterate over every DGSEM element
        @batch for element in eachelement(dg, cache)
            # iterate over every node
            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, element)

                # keep rho_rain positive
                if (u_node[4] < 0.0)
                    u[4, i, j, element] = 0.0
                end

                # keep rho_cloud positive
                if (u_node[3] < 0.0)
                    u[3, i, j, element] = 0.0
                end
            end
        end
    end
end #muladd end
