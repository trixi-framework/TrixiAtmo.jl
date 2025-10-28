@muladd begin
#! format: noindnent

function init_auxiliary_node_variables!(aux_values, mesh::DGMultiMesh,
                                        equations::AbstractCovariantEquations{2, 3},
                                        dg::DGMulti{<:Any, <:Wedge},
                                        bottom_topography)
    rd = dg.basis
    (; V1) = rd
    (; xyz) = mesh.md
    md = mesh.md
    n_aux = n_aux_node_vars(equations)

    # Compute the radius of the sphere from the first element's fourth vertex, such that we can use it
    # throughout the computation. We assume that each Wedge element's last three corner vertices lie
    # on the simulated sphere.
    VX, VY, VZ = map(coords -> transpose(coords[:, 1]) / V1', xyz)
    v_outer = getindex.([VX, VY, VZ], 4)
    radius = norm(v_outer)

    for element in eachelement(mesh, dg)
        # Compute corner vertices of the element
        VX, VY, VZ = map(coords -> transpose(coords[:, element]) / V1', xyz)
        v1, v2, v3 = map(i -> getindex.([VX, VY, VZ], i), 4:6)

        aux_node = Vector{eltype(aux_values[1, 1])}(undef, n_aux)
        
        # Compute the auxiliary metric information at each node
        for i in 1:Trixi.nnodes(dg)
            # Covariant basis in the desired global coordinate system as columns of a matrix
            basis_covariant = calc_basis_covariant(v1, v2, v3,
                                                   rd.rst[1][i], rd.rst[2][i],
                                                   radius,
                                                   equations.global_coordinate_system)
            
            aux_node[1:6] = SVector(basis_covariant)
            
            
            # Covariant metric tensor G := basis_covariant' * basis_covariant
            metric_covariant = basis_covariant' * basis_covariant

            # Contravariant metric tensor inv(G)
            metric_contravariant = inv(metric_covariant)

            # Contravariant basis vectors as rows of a matrix
            basis_contravariant = metric_contravariant * basis_covariant'


            aux_node[7:12] = SVector(basis_contravariant)
            # Area element
            aux_node[13] = sqrt(det(metric_covariant))

            # Covariant metric tensor components
            aux_node[14:16] = SVector(metric_covariant[1, 1],
                                                          metric_covariant[1, 2],
                                                          metric_covariant[2, 2])

            # Contravariant metric tensor components
            aux_node[17:19] = SVector(metric_contravariant[1, 1],
                                                          metric_contravariant[1, 2],
                                                          metric_contravariant[2, 2])
            # Bottom topography
            if !isnothing(bottom_topography)
                nothing
            end
            aux_values[i, element] = SVector{n_aux}(aux_node)
        end
        # TODO: implement Christoffel symbols
        # Christoffel symbols of the second kind (aux_values[21:26, :, :, element])
        # calc_christoffel_symbols!(aux_values, mesh, equations, dg, element)
    end

    return nothing
end

# Analytically compute the transformation matrix A, such that G = AᵀA is the 
# covariant metric tensor and a_i = A[1,i] * e_x + A[2,i] * e_y + A[3,i] * e_z denotes 
# the covariant tangent basis, where e_x, e_y, and e_z are the Cartesian unit basis vectors.
@inline function calc_basis_covariant(v1, v2, v3, xi1, xi2, radius,
                                      ::GlobalCartesianCoordinates)

    # Construct a bilinear mapping based on the four corner vertices
    xe = 0.5f0 * (-(xi1 + xi2) * v1 + (1 + xi1) * v2 +
          (1 + xi2) * v3)
    # Derivatives of bilinear map with respect to reference coordinates xi1, xi2
    dxedxi1 = 0.5f0 *
              (-v1 + v2)
    dxedxi2 = 0.5f0 *
              (-v1 + v3)

    # Use product/quotient rule on the projection
    norm_xe = norm(xe)
    dxdxi1 = radius / norm_xe * (dxedxi1 - dot(xe, dxedxi1) / norm_xe^2 * xe)
    dxdxi2 = radius / norm_xe * (dxedxi2 - dot(xe, dxedxi2) / norm_xe^2 * xe)

    return SMatrix{3, 2}(dxdxi1[1], dxdxi1[2], dxdxi1[3],
                         dxdxi2[1], dxdxi2[2], dxdxi2[3])
end
end # @muladd
