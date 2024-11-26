@muladd begin
#! format: noindent

# Using element_local mapping
function P4estMeshQuadIcosahedron2D(radius;
                                    polydeg, RealT = Float64,
                                    initial_refinement_level = 0,
                                    unsaved_changes = true,
                                    p4est_partition_allow_for_coarsening = true)
    connectivity = connectivity_icosahedron_2D()

    n_trees = 60 # 20 triangles subdivided into 3 quads each

    basis = LobattoLegendreBasis(RealT, polydeg)
    nodes = basis.nodes

    tree_node_coordinates = Array{RealT, 4}(undef, 3,
                                            ntuple(_ -> length(nodes), 2)...,
                                            n_trees)
    calc_tree_node_coordinates_quad_icosahedron_local!(tree_node_coordinates, nodes,
                                                       radius)

    p4est = Trixi.new_p4est(connectivity, initial_refinement_level)

    boundary_names = fill(Symbol("---"), 2 * 2, n_trees)

    return P4estMesh{2}(p4est, tree_node_coordinates, nodes,
                        boundary_names, "", unsaved_changes,
                        p4est_partition_allow_for_coarsening)
end

function connectivity_icosahedron_2D()
    # Illustration of the unfolded icosahedron with the numbering of the triangular faces
    #
    #           ,'|.     
    #         ,'  | `.   
    #       ,' 4  | 3 `. 
    #     ,'_____ |_____`.
    #     \      /\      /
    #      \ 5  /  \ 2  /
    #       \  /  1 \  /
    #        \/_______/_______________________________     
    #         \      /\      /\      /\      /\      /\    
    #          \ 6  /  \ 7  /  \ 8  /  \ 9  /  \ 10 /  \   
    #           \  / 11 \  / 12 \  / 13 \  / 14 \  / 15 \  
    #            \/_______/_______/_______/_______/______\
    #                                            /\      /\   
    #                                           /  \ 20 /  \  
    #                                          / 19 \  / 16 \ 
    #                                         /______\/______\
    #                                          .     |      ,'
    #                                           `.18 | 17 ,'  
    #                                             `. |  ,'    
    #                                               `|,'       
    #
    # Each triangle is subdivided into 3 quadrilaterals with the local coordinate system:
    #
    #                      /\
    #                     /  \
    #                    /    \
    #                   /      \
    #                  /        \
    #                 /          \
    #                /     3      \
    #               /              \
    #              /  η         ξ   \
    #             /   ⎻⎼⎼⎽⎽ ⎽⎽⎼⎼⎻    \
    #            /‾⎺⎺⎻⎻⎼⎼⎽⎽ ⎽⎽⎼⎼⎻⎻⎺⎺‾ \
    #           /          |           \
    #          /           |            \
    #         /            |             \
    #        /   η         | ↑η           \
    #       /  /     1     | |      2      \ 
    #      /  /            | |              \
    #     /  /             | |               \
    #    /  -------->ξ     | └------->ξ       \
    #   /__________________|___________________\

    num_triangles = 20
    trees_per_triangle = 3

    linear_indices = LinearIndices((trees_per_triangle, num_triangles))

    # Vertices represent the coordinates of the forest. This is used by `p4est`
    # to write VTK files.
    # Trixi.jl doesn't use the coordinates from `p4est`, so the vertices can be empty.
    n_vertices = 0
    n_trees = trees_per_triangle * num_triangles

    # No corner connectivity is needed
    n_corners = 0
    vertices = Trixi.C_NULL
    tree_to_vertex = Trixi.C_NULL

    tree_to_tree = Array{Trixi.p4est_topidx_t, 2}(undef, 4, n_trees)
    tree_to_face = Array{Int8, 2}(undef, 4, n_trees)

    # Connectivities for the first 5 (1:5) and the last 5 (16:20) triangular faces
    # Notes:
    # - We subtract 1 because `p4est` uses zero-based indexing)
    # - We use circshift to do a circular shift in the triangle list of the 5 triangles that 
    #   are connected (useful due to periodicity)
    # - We use triangle_offset_base to store the offset for the connections along the base of the 
    #   triangles. For example, triangle 1 is connected to 6 (offset = +5), and triangle 16
    #   is connected to 11 (offset = -5)
    # - We use triangle_offset_13 to store the offset for the connections along the edge of the 
    #   triangles where quads 1 and 3 are, and triangle_offset_23 to store the offset for the connections
    #   along the edge of the triangles where quads 1 and 2 are. For example, triangle 1 is connected to 5 
    #   along the side where quads 1 and 3 are (triangle_offset_13 = +4), and triangle 1 is also 
    #   connected to 2 along the side where quads 2 and 3 are (triangle_offset_23 = 1)
    triangle_list = Vector(101:105)
    triangle_list_1 = Vector(1:5)
    triangle_list_2 = Vector(16:20)
    # initialize triangle_offset_base with a random value (scopes!)
    triangle_offset_base = 1000

    # Loop over the triangles
    for triangle in [1:5; 16:20]
        if triangle <= 5
            triangle_offset_base = 5
            triangle_list = triangle_list_1
            triangle_offset_13 = 4
            triangle_offset_23 = 1
            offset = 0
        else # triangle >=16
            triangle_offset_base = -5
            triangle_list = triangle_list_2
            triangle_offset_13 = 1
            triangle_offset_23 = 4
            offset = 15
        end

        # Quad 1
        ########
        tree = linear_indices[1, triangle]

        # Negative xi-direction
        tree_to_tree[1, tree] = linear_indices[2,
                                               circshift(triangle_list,
                                                         -triangle_offset_13)[triangle - offset]] -
                                1
        tree_to_face[1, tree] = 1

        # Positive xi-direction
        tree_to_tree[2, tree] = linear_indices[2, triangle] - 1
        tree_to_face[2, tree] = 0

        # Negative eta-direction
        tree_to_tree[3, tree] = linear_indices[2, triangle + triangle_offset_base] - 1
        tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4

        # Positive eta-direction
        tree_to_tree[4, tree] = linear_indices[3, triangle] - 1
        tree_to_face[4, tree] = 4 # first face dimensions are oppositely oriented, add 4

        # Quad 2
        ########
        tree = linear_indices[2, triangle]

        # Negative xi-direction
        tree_to_tree[1, tree] = linear_indices[1, triangle] - 1
        tree_to_face[1, tree] = 1

        # Positive xi-direction
        tree_to_tree[2, tree] = linear_indices[1,
                                               circshift(triangle_list,
                                                         -triangle_offset_23)[triangle - offset]] -
                                1
        tree_to_face[2, tree] = 0

        # Negative eta-direction
        tree_to_tree[3, tree] = linear_indices[1, triangle + triangle_offset_base] - 1
        tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4

        # Positive eta-direction
        tree_to_tree[4, tree] = linear_indices[3, triangle] - 1
        tree_to_face[4, tree] = 2

        # Quad 3
        ########
        tree = linear_indices[3, triangle]

        # Negative xi-direction
        tree_to_tree[1, tree] = linear_indices[1, triangle] - 1
        tree_to_face[1, tree] = 7 # first face dimensions are oppositely oriented, add 4

        # Positive xi-direction
        tree_to_tree[2, tree] = linear_indices[3,
                                               circshift(triangle_list,
                                                         -triangle_offset_23)[triangle - offset]] -
                                1
        tree_to_face[2, tree] = 3

        # Negative eta-direction
        tree_to_tree[3, tree] = linear_indices[2, triangle] - 1
        tree_to_face[3, tree] = 3

        # Positive eta-direction
        tree_to_tree[4, tree] = linear_indices[3,
                                               circshift(triangle_list,
                                                         -triangle_offset_13)[triangle - offset]] -
                                1
        tree_to_face[4, tree] = 1
    end

    # Connectivities for the triangular faces 6:15
    # Notes
    # - We subtract 1 because `p4est` uses zero-based indexing)
    # - We use triangle_connectivity_n* to store the triangle connectivity along edges
    # - We use circshift to do a circular shift in the triangle list of triangles 6:15
    #   (useful due to periodicity)
    # - We use triangle_offset_base to store the offset for the connections along the base of the 
    #   triangles. For example, triangle 6 is connected to 1 (offset = -5), and triangle 11
    #   is connected to 16 (offset = +5)
    triangle_list = Vector(6:15)
    # Triangle connectivity along edges with orientation north-east (/)
    triangle_connectivity_ne = circshift(triangle_list, -5)
    # Triangle connectivity along edges with orientation north-west (\)
    triangle_connectivity_nw = [15; 11:14; 7:10; 6]
    # initialize triangle_offset_base (scopes!)
    triangle_offset_base = 1000

    # Loop over the triangles
    for triangle in 6:15
        if triangle <= 10
            triangle_offset_base = -5
        else # triangle >=11
            triangle_offset_base = 5
        end

        # Quad 1
        ########
        tree = linear_indices[1, triangle]

        # Negative xi-direction
        tree_to_tree[1, tree] = linear_indices[3,
                                               triangle_connectivity_ne[triangle - 5]] -
                                1
        tree_to_face[1, tree] = 7 # first face dimensions are oppositely oriented, add 4

        # Positive xi-direction
        tree_to_tree[2, tree] = linear_indices[2, triangle] - 1
        tree_to_face[2, tree] = 0

        # Negative eta-direction
        tree_to_tree[3, tree] = linear_indices[2, triangle + triangle_offset_base] - 1
        tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4

        # Positive eta-direction
        tree_to_tree[4, tree] = linear_indices[3, triangle] - 1
        tree_to_face[4, tree] = 4 # first face dimensions are oppositely oriented, add 4

        # Quad 2
        ########
        tree = linear_indices[2, triangle]

        # Negative xi-direction
        tree_to_tree[1, tree] = linear_indices[1, triangle] - 1
        tree_to_face[1, tree] = 1

        # Positive xi-direction
        tree_to_tree[2, tree] = linear_indices[3,
                                               triangle_connectivity_nw[triangle - 5]] -
                                1
        tree_to_face[2, tree] = 5 # first face dimensions are oppositely oriented, add 4

        # Negative eta-direction
        tree_to_tree[3, tree] = linear_indices[1, triangle + triangle_offset_base] - 1
        tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4

        # Positive eta-direction
        tree_to_tree[4, tree] = linear_indices[3, triangle] - 1
        tree_to_face[4, tree] = 2

        # Quad 3
        ########
        tree = linear_indices[3, triangle]

        # Negative xi-direction
        tree_to_tree[1, tree] = linear_indices[1, triangle] - 1
        tree_to_face[1, tree] = 7 # first face dimensions are oppositely oriented, add 4

        # Positive xi-direction
        tree_to_tree[2, tree] = linear_indices[2,
                                               triangle_connectivity_nw[triangle - 5]] -
                                1
        tree_to_face[2, tree] = 5 # first face dimensions are oppositely oriented, add 4

        # Negative eta-direction
        tree_to_tree[3, tree] = linear_indices[2, triangle] - 1
        tree_to_face[3, tree] = 3

        # Positive eta-direction
        tree_to_tree[4, tree] = linear_indices[1,
                                               triangle_connectivity_ne[triangle - 5]] -
                                1
        tree_to_face[4, tree] = 4 # first face dimensions are oppositely oriented, add 4
    end

    tree_to_corner = Trixi.C_NULL
    # `p4est` docs: "in trivial cases it is just a pointer to a p4est_topix value of 0."
    # We don't need corner connectivity, so this is a trivial case.
    ctt_offset = zeros(Trixi.p4est_topidx_t, 1)

    corner_to_tree = Trixi.C_NULL
    corner_to_corner = Trixi.C_NULL

    connectivity = Trixi.p4est_connectivity_new_copy(n_vertices, n_trees, n_corners,
                                                     vertices, tree_to_vertex,
                                                     tree_to_tree, tree_to_face,
                                                     tree_to_corner, ctt_offset,
                                                     corner_to_tree, corner_to_corner)

    @assert Trixi.p4est_connectivity_is_valid(connectivity) == 1

    return connectivity
end

# Calculate physical coordinates of each node of a 2D quad-icosahedral mesh using the
# element-local mapping from Guba et al. (see https://doi.org/10.5194/gmd-7-2803-2014,
# Appendix A). 
function calc_tree_node_coordinates_quad_icosahedron_local!(node_coordinates::AbstractArray{<:Any,
                                                                                            4},
                                                            nodes, radius)
    num_triangles = 20
    trees_per_triangle = 3

    linear_indices = LinearIndices((trees_per_triangle, num_triangles))

    triangle_vertices = Array{eltype(node_coordinates), 2}(undef, 3, 7)
    icosahedron_vertices = calc_icosahedron_vertices(radius)

    for triangle in 1:num_triangles
        calc_icosahedron_triangle_vertices!(triangle_vertices, icosahedron_vertices,
                                            radius, triangle)

        for local_tree in 1:3
            tree = linear_indices[local_tree, triangle]

            idx = tree_vertices_idx(local_tree)

            # Vertices for bilinear mapping
            v1 = triangle_vertices[:, idx[1]]
            v2 = triangle_vertices[:, idx[2]]
            v3 = triangle_vertices[:, idx[3]]
            v4 = triangle_vertices[:, idx[4]]

            for j in eachindex(nodes), i in eachindex(nodes)
                # node_coordinates are the mapped reference node coordinates
                node_coordinates[:, i, j, tree] .= local_mapping(nodes[i], nodes[j],
                                                                 v1, v2, v3, v4, radius)
            end
        end
    end
end

# Icosahedron vertices
#            5
#           ,'|.     
#         ,'  | `.   
#       ,'    |   `. 
#    6,'_____1|_____`.4
#     \      /\      /
#      \    /  \    /
#       \  /    \  /
#       2\/______3/______4_______5_______6_______2     
#         \      /\      /\      /\      /\      /\    
#          \    /  \    /  \    /  \    /  \    /  \   
#           \  /    \  /    \  /    \  /    \  /    \  
#            \/_______/_______/_______/_______/______\7
#            7       8       9      10      11\      /\   
#                                           /  \    /  \  
#                                          /    \  /    \ 
#                                         /______\/______\
#                                       10 .     |12    ,'8
#                                           `.   |    ,'  
#                                             `. |  ,'    
#                                               `|,'
#                                                9
function calc_icosahedron_vertices(radius; RealT = Float64)
    vertices = Array{RealT, 2}(undef, 3, 12)

    vertices[:, 1] = [0, 0, 1]
    vertices[:, 2] = [
        -sqrt(1 / 10 * (5 - sqrt(5))),
        -1 / 2 - 1 / (2 * sqrt(5)),
        1 / sqrt(5)
    ]
    vertices[:, 3] = [
        sqrt(1 / 10 * (5 - sqrt(5))),
        -1 / 2 - 1 / (2 * sqrt(5)),
        1 / sqrt(5)
    ]
    vertices[:, 4] = [
        sqrt(1 / 10 * (5 + sqrt(5))),
        1 / 2 - 1 / (2 * sqrt(5)),
        1 / sqrt(5)
    ]
    vertices[:, 5] = [0, 2 / sqrt(5), 1 / sqrt(5)]
    vertices[:, 6] = [
        -sqrt(1 / 10 * (5 + sqrt(5))),
        1 / 2 - 1 / (2 * sqrt(5)),
        1 / sqrt(5)
    ]
    vertices[:, 7] = [0, -2 / sqrt(5), -1 / sqrt(5)]
    vertices[:, 8] = [
        sqrt(1 / 10 * (5 + sqrt(5))),
        1 / (2 * sqrt(5)) - 1 / 2,
        -1 / sqrt(5)
    ]
    vertices[:, 9] = [
        sqrt(1 / 10 * (5 - sqrt(5))),
        1 / 2 + 1 / (2 * sqrt(5)),
        -1 / sqrt(5)
    ]
    vertices[:, 10] = [
        -sqrt(1 / 10 * (5 - sqrt(5))),
        1 / 2 + 1 / (2 * sqrt(5)),
        -1 / sqrt(5)
    ]
    vertices[:, 11] = [
        -sqrt(1 / 10 * (5 + sqrt(5))),
        1 / (2 * sqrt(5)) - 1 / 2,
        -1 / sqrt(5)
    ]
    vertices[:, 12] = [0, 0, -1]

    return vertices * radius
end

function icosahedron_triangle_vertices_idx(triangle)
    if triangle == 1
        return [2, 3, 1]
    elseif triangle == 2
        return [3, 4, 1]
    elseif triangle == 3
        return [4, 5, 1]
    elseif triangle == 4
        return [5, 6, 1]
    elseif triangle == 5
        return [6, 2, 1]
    elseif triangle == 6
        return [3, 2, 7]
    elseif triangle == 7
        return [4, 3, 8]
    elseif triangle == 8
        return [5, 4, 9]
    elseif triangle == 9
        return [6, 5, 10]
    elseif triangle == 10
        return [2, 6, 11]
    elseif triangle == 11
        return [7, 8, 3]
    elseif triangle == 12
        return [8, 9, 4]
    elseif triangle == 13
        return [9, 10, 5]
    elseif triangle == 14
        return [10, 11, 6]
    elseif triangle == 15
        return [11, 7, 2]
    elseif triangle == 16
        return [8, 7, 12]
    elseif triangle == 17
        return [9, 8, 12]
    elseif triangle == 18
        return [10, 9, 12]
    elseif triangle == 19
        return [11, 10, 12]
    elseif triangle == 20
        return [7, 11, 12]
    end
end

function calc_icosahedron_triangle_vertices!(triangle_vertices, icosahedron_vertices,
                                             radius, triangle)
    # Retrieve triangle vertives
    corners_triangle = icosahedron_vertices[:,
                                            icosahedron_triangle_vertices_idx(triangle)]

    # local vertex numbering for triangle (corners in parenthesis)
    #                       5 (3)
    #                      /\
    #                     /  \
    #                    /    \
    #                   /      \
    #                  /        \
    #                 /          \
    #                /            \
    #               /              \
    #              /                \
    #            6/                 4\
    #            /⎺⎻⎼⎽          ⎽⎼⎻⎺  \
    #           /     ⎺⎻⎼⎽7 ⎽⎼⎻⎺       \
    #          /          ⎺|            \
    #         /            |             \
    #        /             |              \
    #       /              |               \ 
    #      /               |                \
    #     /                |                 \
    #    /                 |                  \
    #  1/_________________2|__________________3\
    #   (1)                                     (2)
    triangle_vertices[:, 1] = corners_triangle[:, 1]

    v2_bilinear = 0.5 * (corners_triangle[:, 1] + corners_triangle[:, 2])
    triangle_vertices[:, 2] = radius * v2_bilinear / norm(v2_bilinear)

    triangle_vertices[:, 3] = corners_triangle[:, 2]

    v4_bilinear = 0.5 * (corners_triangle[:, 2] + corners_triangle[:, 3])
    triangle_vertices[:, 4] = radius * v4_bilinear / norm(v4_bilinear)

    triangle_vertices[:, 5] = corners_triangle[:, 3]

    v6_bilinear = 0.5 * (corners_triangle[:, 1] + corners_triangle[:, 3])
    triangle_vertices[:, 6] = radius * v6_bilinear / norm(v6_bilinear)

    v7_bilinear = (corners_triangle[:, 1] + corners_triangle[:, 2] +
                   corners_triangle[:, 3]) / 3
    triangle_vertices[:, 7] = radius * v7_bilinear / norm(v7_bilinear)
end

function tree_vertices_idx(local_tree)
    if local_tree == 1
        return [1, 2, 7, 6]
    elseif local_tree == 2
        return [2, 3, 4, 7]
    else
        return [7, 4, 5, 6]
    end
end
end
