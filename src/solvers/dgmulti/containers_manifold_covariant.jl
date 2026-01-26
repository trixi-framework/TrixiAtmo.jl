@muladd begin
#! format: noindnent

function init_auxiliary_node_variables!(aux_values, aux_quad_values, aux_face_values, mesh::DGMultiMesh,
                                        equations::AbstractCovariantEquations{2, 3},
                                        dg::DGMulti{<:Any, <:Tri},
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
    v_outer = getindex.([VX, VY, VZ], 1)
    radius = norm(v_outer)

    for element in eachelement(mesh, dg)
        # Compute corner vertices of the element
        VX, VY, VZ = map(coords -> transpose(coords[:, element]) / V1', xyz)
        v1, v2, v3 = map(i -> getindex.([VX, VY, VZ], i), 1:3)

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
                x_node = map(coords -> coords[i, element], xyz)
                aux_node[20] = bottom_topography(x_node)
            else
                aux_node[20] = zero(eltype(aux_node))
            end
            aux_values[i, element] = SVector{n_aux}(aux_node)
        end
        # Christoffel symbols of the second kind (aux_values[21:26, :, :, element])
        calc_christoffel_symbols!(aux_values, mesh, equations, dg, element)

        aux_node = Vector{eltype(aux_values[1, 1])}(undef, n_aux)
        
        # Compute the auxiliary metric information at each quad node
        for i in Trixi.each_quad_node(mesh, dg)
            # Covariant basis in the desired global coordinate system as columns of a matrix
            basis_covariant = calc_basis_covariant(v1, v2, v3,
                                                   rd.rstq[1][i], rd.rstq[2][i],
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
                x_node = map(coords -> coords[i, element], xyz)
                aux_node[20] = bottom_topography(x_node)
            else
                aux_node[20] = zero(eltype(aux_node))
            end
            aux_quad_values[i, element] = SVector{n_aux}(aux_node)
        end

        aux_node = Vector{eltype(aux_values[1, 1])}(undef, n_aux)
        
        # Compute the auxiliary metric information at each face node
        for i in Trixi.each_face_node(mesh, dg)
            # Covariant basis in the desired global coordinate system as columns of a matrix
            basis_covariant = calc_basis_covariant(v1, v2, v3,
                                                   rd.rstf[1][i], rd.rstf[2][i],
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
                x_node = map(coords -> coords[i, element], xyz)
                aux_node[20] = bottom_topography(x_node)
            else
                aux_node[20] = zero(eltype(aux_node))
            end
            aux_face_values[i, element] = SVector{n_aux}(aux_node)
        end
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

function calc_christoffel_symbols!(aux_values, mesh::DGMultiMesh,
                                   equations::AbstractCovariantEquations{2, 3}, dg,
                                   element)
    rd = dg.basis
    (; Vq, Drst) = rd
    Dr, Ds = Drst

    for i in 1:Trixi.nnodes(dg)

        # Numerically differentiate covariant metric components with respect to ξ¹
        dG11dxi1 = zero(eltype(aux_values[i, element]))
        dG12dxi1 = zero(eltype(aux_values[i, element]))
        dG22dxi1 = zero(eltype(aux_values[i, element]))
        for jj in 1:nnodes(dg)
            aux_node_jj = aux_values[jj, element]
            Gcov_jj = metric_covariant(aux_node_jj, equations)
            dG11dxi1 += Dr[i, jj] * Gcov_jj[1, 1]
            dG12dxi1 += Dr[i, jj] * Gcov_jj[1, 2]
            dG22dxi1 += Dr[i, jj] * Gcov_jj[2, 2]
        end

        # Numerically differentiate covariant metric components with respect to ξ²
        dG11dxi2 = zero(eltype(aux_values[i, element]))
        dG12dxi2 = zero(eltype(aux_values[i, element]))
        dG22dxi2 = zero(eltype(aux_values[i, element]))
        for jj in 1:nnodes(dg)
            aux_node_jj = aux_values[jj, element]
            Gcov_jj = metric_covariant(aux_node_jj, equations)
            dG11dxi2 += Ds[i, jj] * Gcov_jj[1, 1]
            dG12dxi2 += Ds[i, jj] * Gcov_jj[1, 2]
            dG22dxi2 += Ds[i, jj] * Gcov_jj[2, 2]
        end

        # Compute Christoffel symbols of the first kind
        christoffel_firstkind_1 = SMatrix{2, 2}(0.5f0 * dG11dxi1,
                                                0.5f0 * dG11dxi2,
                                                0.5f0 * dG11dxi2,
                                                dG12dxi2 - 0.5f0 * dG22dxi1)
        christoffel_firstkind_2 = SMatrix{2, 2}(dG12dxi1 - 0.5f0 * dG11dxi2,
                                                0.5f0 * dG22dxi1,
                                                0.5f0 * dG22dxi1,
                                                0.5f0 * dG22dxi2)

        # Raise indices to get Christoffel symbols of the second kind
        aux_node = Vector(aux_values[i, element])
        Gcon = metric_contravariant(aux_node, equations)
        aux_node[21] = Gcon[1, 1] * christoffel_firstkind_1[1, 1] +
                                           Gcon[1, 2] * christoffel_firstkind_2[1, 1]
        aux_node[22] = Gcon[1, 1] * christoffel_firstkind_1[1, 2] +
                                           Gcon[1, 2] * christoffel_firstkind_2[1, 2]
        aux_node[23] = Gcon[1, 1] * christoffel_firstkind_1[2, 2] +
                                           Gcon[1, 2] * christoffel_firstkind_2[2, 2]

        aux_node[24] = Gcon[2, 1] * christoffel_firstkind_1[1, 1] +
                                           Gcon[2, 2] * christoffel_firstkind_2[1, 1]
        aux_node[25] = Gcon[2, 1] * christoffel_firstkind_1[1, 2] +
                                           Gcon[2, 2] * christoffel_firstkind_2[1, 2]
        aux_node[26] = Gcon[2, 1] * christoffel_firstkind_1[2, 2] +
                                           Gcon[2, 2] * christoffel_firstkind_2[2, 2]

        aux_values[i, element] = SVector{n_aux_node_vars(equations)}(aux_node)
    end
end
end # @muladd
