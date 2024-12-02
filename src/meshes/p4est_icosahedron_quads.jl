@muladd begin
#! format: noindent

# Using element_local mapping
"""
    P4estMeshQuadIcosahedron2D(trees_per_face_dimension, radius;
                               polydeg, RealT=Float64,
                               initial_refinement_level=0, unsaved_changes=true,
                               p4est_partition_allow_for_coarsening=true)

Build a quad-based icosahedral mesh as a 2D `P4estMesh` with
`60 * trees_per_face_dimension^2` trees (20 triangular faces of the icosahedron,
each subdivided into 3 parent quads, each of which subdivided into `trees_per_face_dimension^2` trees).

The node coordinates of the trees will be obtained using the element-local mapping from
Appendix A of [Guba et al. (2014)](https://doi.org/10.5194/gmd-7-2803-2014). 
See [`P4estMeshCubedSphere2D`](@ref) for more information about the element-local mapping.

The mesh will have no boundaries.

# Arguments
- `trees_per_face_dimension::Integer`: the number of trees in the two local dimensions of
                                       each parent quad.
- `radius::Integer`: the radius of the sphere.
- `polydeg::Integer`: polynomial degree used to store the geometry of the mesh.
                      The mapping will be approximated by an interpolation polynomial
                      of the specified degree for each tree.
- `RealT::Type`: the type that should be used for coordinates.
- `initial_refinement_level::Integer`: refine the mesh uniformly to this level before the
  simulation starts.
- `unsaved_changes::Bool`: if set to `true`, the mesh will be saved to a mesh file.
- `p4est_partition_allow_for_coarsening::Bool`: Must be `true` when using AMR to make mesh
  adaptivity independent of domain partitioning. Should be `false` for static meshes to
  permit more fine-grained partitioning.

!!! warning 
    Adaptivity and MPI parallelization are not yet supported for equations in covariant 
    form, and we require `initial_refinement_level = 0` for such cases. Furthermore, the 
    calculation of the metric terms for the covariant form currently requires `polydeg` to 
    be equal to the polynomial degree of the solver.
!!!
"""
function P4estMeshQuadIcosahedron2D(trees_per_face_dimension, radius;
                                    polydeg, RealT = Float64,
                                    initial_refinement_level = 0,
                                    unsaved_changes = true,
                                    p4est_partition_allow_for_coarsening = true)
    connectivity = connectivity_icosahedron_2D(trees_per_face_dimension)

    n_trees = 60 * trees_per_face_dimension^2 # 20 triangles subdivided into 3 quads each

    basis = LobattoLegendreBasis(RealT, polydeg)
    nodes = basis.nodes

    tree_node_coordinates = Array{RealT, 4}(undef, 3,
                                            ntuple(_ -> length(nodes), 2)...,
                                            n_trees)
    calc_tree_node_coordinates_quad_icosahedron_local!(tree_node_coordinates, nodes,
                                                       trees_per_face_dimension, radius)

    p4est = Trixi.new_p4est(connectivity, initial_refinement_level)

    boundary_names = fill(Symbol("---"), 2 * 2, n_trees)

    return P4estMesh{2}(p4est, tree_node_coordinates, nodes,
                        boundary_names, "", unsaved_changes,
                        p4est_partition_allow_for_coarsening)
end

# Fig 1: Illustration of the unfolded icosahedron with the numbering of the triangular faces
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
# Each triangle is subdivided into 3 quadrilaterals with a local (ξ,η)-coordinate system.
#
# Fig 2:
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
#
# Each of those quadrilaterlas is subdivided into trees_per_face_dimension^2 trees.
#
# We use the following numbering for the 12 vertices of the icosahedron
# Fig 3:
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

# Function to compute the vertices' coordinates of an icosahedron inscribed in a sphere of radius `radius`
function calc_node_coordinates_icosahedron_vertices(radius; RealT = Float64)
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

# Index map for the corner vertices of the triangular faces on the icosahedron (see Fig 1 and Fig 3)
# We use a counter-clockwise numbering
const icosahedron_triangle_vertices_idx_map = ([2, 3, 1], # Triangle 1
                                               [3, 4, 1], # Triangle 2
                                               [4, 5, 1], # Triangle 3
                                               [5, 6, 1], # Triangle 4
                                               [6, 2, 1], # Triangle 5
                                               [3, 2, 7], # Triangle 6
                                               [4, 3, 8], # Triangle 7
                                               [5, 4, 9], # Triangle 8
                                               [6, 5, 10], # Triangle 9
                                               [2, 6, 11], # Triangle 10
                                               [7, 8, 3], # Triangle 11
                                               [8, 9, 4], # Triangle 12
                                               [9, 10, 5], # Triangle 13
                                               [10, 11, 6], # Triangle 14
                                               [11, 7, 2], # Triangle 15
                                               [8, 7, 12], # Triangle 16
                                               [9, 8, 12], # Triangle 17
                                               [10, 9, 12], # Triangle 18
                                               [11, 10, 12], # Triangle 19
                                               [7, 11, 12])

# We use a local numbering to obtain the quad vertices of each triangular face
#
# Fig 4: Local quad vertex numbering for a triangular face (corner vertices of the triangular face in parenthesis)
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

# Index map for the vertices of each quad on the triangular faces of the icosahedron (see Fig 4)
const icosahedron_quad_vertices_idx_map = ([1, 2, 7, 6], # Quad 1
                                           [2, 3, 4, 7], # Quad 2
                                           [7, 4, 5, 6]) # Quad 3
end

# Function to initialize the p4est connectivity for the icosahedral grid. 
# For reference, see Fig 1 and Fig 2 above.
function connectivity_icosahedron_2D(trees_per_face_dimension)
    num_triangles = 20
    trees_per_triangle = 3
    n_cells = trees_per_face_dimension

    linear_indices = LinearIndices((n_cells, n_cells, trees_per_triangle,
                                    num_triangles))

    # Vertices represent the coordinates of the forest. This is used by `p4est`
    # to write VTK files.
    # Trixi.jl doesn't use the coordinates from `p4est`, so the vertices can be empty.
    n_vertices = 0
    n_trees = n_cells * n_cells * trees_per_triangle * num_triangles

    # No corner connectivity is needed
    n_corners = 0
    vertices = Trixi.C_NULL
    tree_to_vertex = Trixi.C_NULL

    tree_to_tree = Array{Trixi.p4est_topidx_t, 2}(undef, 4, n_trees)
    tree_to_face = Array{Int8, 2}(undef, 4, n_trees)

    # Connectivities for the first 5 (1:5) and the last 5 (16:20) triangular faces
    # Notes:
    # - We subtract 1 because `p4est` uses zero-based indexing
    # - We use circshift to do a circular shift in the triangle list of the 5 triangles that 
    #   are connected to each pole (useful due to periodicity)
    # - We use triangle_offset_base to store the offset for the connections along the base of the 
    #   triangles. For example, triangle 1 is connected to 6 (offset = +5), and triangle 16
    #   is connected to 11 (offset = -5)
    # - We use triangle_offset_13 to store the offset for the connections along the edge of the 
    #   triangles where quads 1 and 3 are, and triangle_offset_23 to store the offset for the connections
    #   along the edge of the triangles where quads 1 and 2 are. For example, triangle 1 is connected to 
    #   triangle 5 along the side where quads 1 and 3 are (triangle_offset_13 = +4), and triangle 1 is also 
    #   connected to triangle 2 along the side where quads 2 and 3 are (triangle_offset_23 = 1)
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
        for cell_y in 1:n_cells, cell_x in 1:n_cells
            tree = linear_indices[cell_x, cell_y, 1, triangle]

            # Negative xi-direction
            if cell_x > 1 # Connect to tree at the same quad
                tree_to_tree[1, tree] = linear_indices[cell_x - 1, cell_y, 1,
                                                       triangle] - 1
                tree_to_face[1, tree] = 1
            else
                tree_to_tree[1, tree] = linear_indices[end, cell_y, 2,
                                                       circshift(triangle_list,
                                                                 -triangle_offset_13)[triangle - offset]] -
                                        1
                tree_to_face[1, tree] = 1
            end

            # Positive xi-direction
            if cell_x < n_cells # Connect to tree at the same quad
                tree_to_tree[2, tree] = linear_indices[cell_x + 1, cell_y, 1,
                                                       triangle] - 1
                tree_to_face[2, tree] = 0
            else
                tree_to_tree[2, tree] = linear_indices[1, cell_y, 2, triangle] - 1
                tree_to_face[2, tree] = 0
            end

            # Negative eta-direction
            if cell_y > 1 # Connect to tree at the same quad
                tree_to_tree[3, tree] = linear_indices[cell_x, cell_y - 1, 1,
                                                       triangle] - 1
                tree_to_face[3, tree] = 3
            else
                tree_to_tree[3, tree] = linear_indices[n_cells - cell_x + 1, 1, 2,
                                                       triangle + triangle_offset_base] -
                                        1
                tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4
            end

            # Positive eta-direction
            if cell_y < n_cells # Connect to tree at the same quad
                tree_to_tree[4, tree] = linear_indices[cell_x, cell_y + 1, 1,
                                                       triangle] - 1
                tree_to_face[4, tree] = 2
            else
                tree_to_tree[4, tree] = linear_indices[1, n_cells - cell_x + 1, 3,
                                                       triangle] - 1
                tree_to_face[4, tree] = 4 # first face dimensions are oppositely oriented, add 4
            end
        end

        # Quad 2
        ########
        for cell_y in 1:n_cells, cell_x in 1:n_cells
            tree = linear_indices[cell_x, cell_y, 2, triangle]

            # Negative xi-direction
            if cell_x > 1 # Connect to tree at the same quad
                tree_to_tree[1, tree] = linear_indices[cell_x - 1, cell_y, 2,
                                                       triangle] - 1
                tree_to_face[1, tree] = 1
            else
                tree_to_tree[1, tree] = linear_indices[end, cell_y, 1, triangle] - 1
                tree_to_face[1, tree] = 1
            end

            # Positive xi-direction
            if cell_x < n_cells # Connect to tree at the same quad
                tree_to_tree[2, tree] = linear_indices[cell_x + 1, cell_y, 2,
                                                       triangle] - 1
                tree_to_face[2, tree] = 0
            else
                tree_to_tree[2, tree] = linear_indices[1, cell_y, 1,
                                                       circshift(triangle_list,
                                                                 -triangle_offset_23)[triangle - offset]] -
                                        1
                tree_to_face[2, tree] = 0
            end

            # Negative eta-direction
            if cell_y > 1 # Connect to tree at the same quad
                tree_to_tree[3, tree] = linear_indices[cell_x, cell_y - 1, 2,
                                                       triangle] - 1
                tree_to_face[3, tree] = 3
            else
                tree_to_tree[3, tree] = linear_indices[n_cells - cell_x + 1, 1, 1,
                                                       triangle + triangle_offset_base] -
                                        1
                tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4
            end

            # Positive eta-direction
            if cell_y < n_cells # Connect to tree at the same quad
                tree_to_tree[4, tree] = linear_indices[cell_x, cell_y + 1, 2,
                                                       triangle] - 1
                tree_to_face[4, tree] = 2
            else
                tree_to_tree[4, tree] = linear_indices[cell_x, 1, 3, triangle] - 1
                tree_to_face[4, tree] = 2
            end
        end

        # Quad 3
        ########
        for cell_y in 1:n_cells, cell_x in 1:n_cells
            tree = linear_indices[cell_x, cell_y, 3, triangle]

            # Negative xi-direction
            if cell_x > 1 # Connect to tree at the same quad
                tree_to_tree[1, tree] = linear_indices[cell_x - 1, cell_y, 3,
                                                       triangle] - 1
                tree_to_face[1, tree] = 1
            else
                tree_to_tree[1, tree] = linear_indices[n_cells - cell_y + 1, end, 1,
                                                       triangle] - 1
                tree_to_face[1, tree] = 7 # first face dimensions are oppositely oriented, add 4
            end

            # Positive xi-direction
            if cell_x < n_cells # Connect to tree at the same quad
                tree_to_tree[2, tree] = linear_indices[cell_x + 1, cell_y, 3,
                                                       triangle] - 1
                tree_to_face[2, tree] = 0
            else
                tree_to_tree[2, tree] = linear_indices[cell_y, end, 3,
                                                       circshift(triangle_list,
                                                                 -triangle_offset_23)[triangle - offset]] -
                                        1
                tree_to_face[2, tree] = 3
            end

            # Negative eta-direction
            if cell_y > 1 # Connect to tree at the same quad
                tree_to_tree[3, tree] = linear_indices[cell_x, cell_y - 1, 3,
                                                       triangle] - 1
                tree_to_face[3, tree] = 3
            else
                tree_to_tree[3, tree] = linear_indices[cell_x, end, 2, triangle] - 1
                tree_to_face[3, tree] = 3
            end

            # Positive eta-direction
            if cell_y < n_cells # Connect to tree at the same quad
                tree_to_tree[4, tree] = linear_indices[cell_x, cell_y + 1, 3,
                                                       triangle] - 1
                tree_to_face[4, tree] = 2
            else
                tree_to_tree[4, tree] = linear_indices[end, cell_x, 3,
                                                       circshift(triangle_list,
                                                                 -triangle_offset_13)[triangle - offset]] -
                                        1
                tree_to_face[4, tree] = 1
            end
        end
    end

    # Connectivities for the triangular faces 6:15
    # Notes
    # - We subtract 1 because `p4est` uses zero-based indexing
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
        for cell_y in 1:n_cells, cell_x in 1:n_cells
            tree = linear_indices[cell_x, cell_y, 1, triangle]

            # Negative xi-direction
            if cell_x > 1 # Connect to tree at the same quad
                tree_to_tree[1, tree] = linear_indices[cell_x - 1, cell_y, 1,
                                                       triangle] - 1
                tree_to_face[1, tree] = 1
            else
                tree_to_tree[1, tree] = linear_indices[n_cells - cell_y + 1, end, 3,
                                                       triangle_connectivity_ne[triangle - 5]] -
                                        1
                tree_to_face[1, tree] = 7 # first face dimensions are oppositely oriented, add 4
            end

            # Positive xi-direction
            if cell_x < n_cells # Connect to tree at the same quad
                tree_to_tree[2, tree] = linear_indices[cell_x + 1, cell_y, 1,
                                                       triangle] - 1
                tree_to_face[2, tree] = 0
            else
                tree_to_tree[2, tree] = linear_indices[1, cell_y, 2, triangle] - 1
                tree_to_face[2, tree] = 0
            end

            # Negative eta-direction
            if cell_y > 1 # Connect to tree at the same quad
                tree_to_tree[3, tree] = linear_indices[cell_x, cell_y - 1, 1,
                                                       triangle] - 1
                tree_to_face[3, tree] = 3
            else
                tree_to_tree[3, tree] = linear_indices[n_cells - cell_x + 1, 1, 2,
                                                       triangle + triangle_offset_base] -
                                        1
                tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4
            end

            # Positive eta-direction
            if cell_y < n_cells # Connect to tree at the same quad
                tree_to_tree[4, tree] = linear_indices[cell_x, cell_y + 1, 1,
                                                       triangle] - 1
                tree_to_face[4, tree] = 2
            else
                tree_to_tree[4, tree] = linear_indices[1, n_cells - cell_x + 1, 3,
                                                       triangle] - 1
                tree_to_face[4, tree] = 4 # first face dimensions are oppositely oriented, add 4
            end
        end

        # Quad 2
        ########
        for cell_y in 1:n_cells, cell_x in 1:n_cells
            tree = linear_indices[cell_x, cell_y, 2, triangle]

            # Negative xi-direction
            if cell_x > 1 # Connect to tree at the same quad
                tree_to_tree[1, tree] = linear_indices[cell_x - 1, cell_y, 2,
                                                       triangle] - 1
                tree_to_face[1, tree] = 1
            else
                tree_to_tree[1, tree] = linear_indices[end, cell_y, 1, triangle] - 1
                tree_to_face[1, tree] = 1
            end

            # Positive xi-direction
            if cell_x < n_cells # Connect to tree at the same quad
                tree_to_tree[2, tree] = linear_indices[cell_x + 1, cell_y, 2,
                                                       triangle] - 1
                tree_to_face[2, tree] = 0
            else
                tree_to_tree[2, tree] = linear_indices[end, n_cells - cell_y + 1, 3,
                                                       triangle_connectivity_nw[triangle - 5]] -
                                        1
                tree_to_face[2, tree] = 5 # first face dimensions are oppositely oriented, add 4
            end

            # Negative eta-direction
            if cell_y > 1 # Connect to tree at the same quad
                tree_to_tree[3, tree] = linear_indices[cell_x, cell_y - 1, 2,
                                                       triangle] - 1
                tree_to_face[3, tree] = 3
            else
                tree_to_tree[3, tree] = linear_indices[n_cells - cell_x + 1, 1, 1,
                                                       triangle + triangle_offset_base] -
                                        1
                tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4
            end

            # Positive eta-direction
            if cell_y < n_cells # Connect to tree at the same quad
                tree_to_tree[4, tree] = linear_indices[cell_x, cell_y + 1, 2,
                                                       triangle] - 1
                tree_to_face[4, tree] = 2
            else
                tree_to_tree[4, tree] = linear_indices[cell_x, 1, 3, triangle] - 1
                tree_to_face[4, tree] = 2
            end
        end

        # Quad 3
        ########
        for cell_y in 1:n_cells, cell_x in 1:n_cells
            tree = linear_indices[cell_x, cell_y, 3, triangle]

            # Negative xi-direction
            if cell_x > 1 # Connect to tree at the same quad
                tree_to_tree[1, tree] = linear_indices[cell_x - 1, cell_y, 3,
                                                       triangle] - 1
                tree_to_face[1, tree] = 1
            else
                tree_to_tree[1, tree] = linear_indices[n_cells - cell_y + 1, end, 1,
                                                       triangle] - 1
                tree_to_face[1, tree] = 7 # first face dimensions are oppositely oriented, add 4
            end

            # Positive xi-direction
            if cell_x < n_cells # Connect to tree at the same quad
                tree_to_tree[2, tree] = linear_indices[cell_x + 1, cell_y, 3,
                                                       triangle] - 1
                tree_to_face[2, tree] = 0
            else
                tree_to_tree[2, tree] = linear_indices[end, n_cells - cell_y + 1, 2,
                                                       triangle_connectivity_nw[triangle - 5]] -
                                        1
                tree_to_face[2, tree] = 5 # first face dimensions are oppositely oriented, add 4
            end

            # Negative eta-direction
            if cell_y > 1 # Connect to tree at the same quad
                tree_to_tree[3, tree] = linear_indices[cell_x, cell_y - 1, 3,
                                                       triangle] - 1
                tree_to_face[3, tree] = 3
            else
                tree_to_tree[3, tree] = linear_indices[cell_x, end, 2, triangle] - 1
                tree_to_face[3, tree] = 3
            end

            # Positive eta-direction
            if cell_y < n_cells # Connect to tree at the same quad
                tree_to_tree[4, tree] = linear_indices[cell_x, cell_y + 1, 3,
                                                       triangle] - 1
                tree_to_face[4, tree] = 2
            else
                tree_to_tree[4, tree] = linear_indices[1, n_cells - cell_x + 1, 1,
                                                       triangle_connectivity_ne[triangle - 5]] -
                                        1
                tree_to_face[4, tree] = 4 # first face dimensions are oppositely oriented, add 4
            end
        end
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
                                                            nodes,
                                                            trees_per_face_dimension,
                                                            radius)
    num_triangles = 20
    trees_per_triangle = 3
    n_cells = trees_per_face_dimension

    linear_indices = LinearIndices((n_cells, n_cells, trees_per_triangle,
                                    num_triangles))

    triangle_vertices = Array{eltype(node_coordinates), 2}(undef, 3, 7)
    icosahedron_vertices = calc_node_coordinates_icosahedron_vertices(radius)

    # Get cell length in reference mesh
    dx = 2 / n_cells
    dy = 2 / n_cells

    # Loop over all the triangles of the mesh
    for triangle in 1:num_triangles
        calc_node_coordinates_triangle_vertices!(triangle_vertices,
                                                 icosahedron_vertices,
                                                 radius, triangle)

        # Loop over each parent quad in each triangle
        for local_tree in 1:3
            idx = icosahedron_quad_vertices_idx_map[local_tree]

            # Vertices of the parent quad
            v1_quad = triangle_vertices[:, idx[1]]
            v2_quad = triangle_vertices[:, idx[2]]
            v3_quad = triangle_vertices[:, idx[3]]
            v4_quad = triangle_vertices[:, idx[4]]

            # Loop over the cells/trees in each parent quad
            for cell_y in 1:n_cells, cell_x in 1:n_cells
                tree = linear_indices[cell_x, cell_y, local_tree, triangle]

                x_0 = -1 + (cell_x - 1) * dx
                y_0 = -1 + (cell_y - 1) * dy
                x_1 = -1 + cell_x * dx
                y_1 = -1 + cell_y * dy

                # Obtain the coordinates of the corner nodes for the tree
                v1 = local_mapping(x_0, y_0, v1_quad, v2_quad, v3_quad, v4_quad, radius)
                v2 = local_mapping(x_1, y_0, v1_quad, v2_quad, v3_quad, v4_quad, radius)
                v3 = local_mapping(x_1, y_1, v1_quad, v2_quad, v3_quad, v4_quad, radius)
                v4 = local_mapping(x_0, y_1, v1_quad, v2_quad, v3_quad, v4_quad, radius)

                for j in eachindex(nodes), i in eachindex(nodes)
                    # node_coordinates are the mapped reference node coordinates
                    node_coordinates[:, i, j, tree] .= local_mapping(nodes[i], nodes[j],
                                                                     v1, v2, v3, v4,
                                                                     radius)
                end
            end
        end
    end
end

function calc_node_coordinates_triangle_vertices!(triangle_vertices,
                                                  icosahedron_vertices,
                                                  radius, triangle)
    # Retrieve triangle vertices
    corners_triangle = icosahedron_vertices[:,
                                            icosahedron_triangle_vertices_idx_map[triangle]]

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
