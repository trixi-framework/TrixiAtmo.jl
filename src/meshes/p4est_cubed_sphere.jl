@muladd begin
#! format: noindent

"""
    P4estMeshCubedSphere2D(trees_per_face_dimension, radius;
                            polydeg, RealT=Float64,
                            initial_refinement_level=0, unsaved_changes=true,
                            p4est_partition_allow_for_coarsening=true,
                            element_local_mapping=false)

Build a "Cubed Sphere" mesh as a 2D `P4estMesh` with
`6 * trees_per_face_dimension^2` trees.

The mesh will have no boundaries.

# Arguments
- `trees_per_face_dimension::Integer`: the number of trees in the two local dimensions of
                                       each face.
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
- `element_local_mapping::Bool`: option to use the alternative element-local mapping from
  Appendix A of [Guba et al. (2014)](https://doi.org/10.5194/gmd-7-2803-2014). If set to 
  `true`, the four corner vertex positions for each element will be obtained through an 
  equiangular gnomonic projection [(Ronchi et al. 1996)](https://doi.org/10.1006/jcph.1996.
  0047), and the tree node coordinates within the element (i.e. the field 
  `tree_node_coordinates`) will be obtained by first using a bilinear mapping based on the 
  four corner vertices, and then projecting the bilinearly mapped nodes onto the spherical 
  surface by normalizing the resulting Cartesian coordinates and scaling by  `radius`. If 
  set to `false`, the equiangular gnomonic projection will be used for all tree node 
  coordinates.

!!! warning 
    Adaptivity and MPI parallelization are not yet supported for equations in covariant 
    form, and we require `initial_refinement_level = 0` for such cases. Furthermore, the 
    calculation of the metric terms for the covariant form currently requires `polydeg` to 
    be equal to the polynomial degree of the solver, and `element_local_mapping = true`.
"""
function P4estMeshCubedSphere2D(trees_per_face_dimension, radius;
                                polydeg, RealT = Float64,
                                initial_refinement_level = 0,
                                unsaved_changes = true,
                                p4est_partition_allow_for_coarsening = true,
                                element_local_mapping = false)
    connectivity = connectivity_cubed_sphere_2D(trees_per_face_dimension)

    n_trees = 6 * trees_per_face_dimension^2

    basis = LobattoLegendreBasis(RealT, polydeg)
    nodes = basis.nodes

    tree_node_coordinates = Array{RealT, 4}(undef, 3,
                                            ntuple(_ -> length(nodes), 2)...,
                                            n_trees)
    if element_local_mapping
        calc_tree_node_coordinates_cubed_sphere_local!(tree_node_coordinates, nodes,
                                                       trees_per_face_dimension, radius)
    else
        calc_tree_node_coordinates_cubed_sphere_standard!(tree_node_coordinates, nodes,
                                                          trees_per_face_dimension,
                                                          radius)
    end

    p4est = Trixi.new_p4est(connectivity, initial_refinement_level)

    boundary_names = fill(Symbol("---"), 2 * 2, n_trees)

    return P4estMesh{2}(p4est, tree_node_coordinates, nodes,
                        boundary_names, "", unsaved_changes,
                        p4est_partition_allow_for_coarsening)
end

# Connectivity of a 2D cubed sphere
function connectivity_cubed_sphere_2D(trees_per_face_dimension)
    n_cells_x = n_cells_y = trees_per_face_dimension

    linear_indices = LinearIndices((trees_per_face_dimension, trees_per_face_dimension,
                                    6))

    # Vertices represent the coordinates of the forest. This is used by `p4est`
    # to write VTK files.
    # Trixi.jl doesn't use the coordinates from `p4est`, so the vertices can be empty.
    n_vertices = 0
    n_trees = 6 * n_cells_x * n_cells_y

    # No corner connectivity is needed
    n_corners = 0
    vertices = Trixi.C_NULL
    tree_to_vertex = Trixi.C_NULL

    tree_to_tree = Array{Trixi.p4est_topidx_t, 2}(undef, 4, n_trees)
    tree_to_face = Array{Int8, 2}(undef, 4, n_trees)

    # Illustration of the local coordinates of each face. ξ and η are the first
    # local coordinates of each face. The third local coordinate ζ is always
    # pointing outwards, which yields a right-handed coordinate system for each face.
    #               ┌────────────────────────────────────────────────────┐
    #              ╱│                                                   ╱│
    #             ╱ │                       ξ <───┐                    ╱ │
    #            ╱  │                            ╱                    ╱  │
    #           ╱   │                4 (+y)     V                    ╱   │
    #          ╱    │                          η                    ╱    │
    #         ╱     │                                              ╱     │
    #        ╱      │                                             ╱      │
    #       ╱       │                                            ╱       │
    #      ╱        │                                           ╱        │
    #     ╱         │                    5 (-z)   η            ╱         │
    #    ╱          │                             ↑           ╱          │
    #   ╱           │                             │          ╱           │
    #  ╱            │                       ξ <───┘         ╱            │
    # ┌────────────────────────────────────────────────────┐    2 (+x)   │
    # │             │                                      │             │
    # │             │                                      │      ξ      │
    # │             │                                      │      ↑      │
    # │    1 (-x)   │                                      │      │      │
    # │             │                                      │      │      │
    # │     ╱│      │                                      │     ╱       │
    # │    V │      │                                      │    V        │
    # │   η  ↓      │                                      │   η         │
    # │      ξ      └──────────────────────────────────────│─────────────┘
    # │            ╱         η   6 (+z)                    │            ╱
    # │           ╱          ↑                             │           ╱
    # │          ╱           │                             │          ╱
    # │         ╱            └───> ξ                       │         ╱
    # │        ╱                                           │        ╱
    # │       ╱                                            │       ╱ Global coordinates:
    # │      ╱                                             │      ╱        y
    # │     ╱                      ┌───> ξ                 │     ╱         ↑
    # │    ╱                      ╱                        │    ╱          │
    # │   ╱                      V      3 (-y)             │   ╱           │
    # │  ╱                      η                          │  ╱            └─────> x
    # │ ╱                                                  │ ╱            ╱
    # │╱                                                   │╱            V
    # └────────────────────────────────────────────────────┘            z
    for direction in 1:6
        for cell_y in 1:n_cells_y, cell_x in 1:n_cells_x
            tree = linear_indices[cell_x, cell_y, direction]

            # Subtract 1 because `p4est` uses zero-based indexing
            # Negative x-direction
            if cell_x > 1 # Connect to tree at the same face
                tree_to_tree[1, tree] = linear_indices[cell_x - 1, cell_y,
                                                       direction] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 1 # This is the -x face
                target = 4
                tree_to_tree[1, tree] = linear_indices[end, cell_y, target] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 2 # This is the +x face
                target = 3
                tree_to_tree[1, tree] = linear_indices[end, cell_y, target] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 3 # This is the -y face
                target = 1
                tree_to_tree[1, tree] = linear_indices[end, cell_y, target] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 4 # This is the +y face
                target = 2
                tree_to_tree[1, tree] = linear_indices[end, cell_y, target] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 5 # This is the -z face
                target = 2
                tree_to_tree[1, tree] = linear_indices[cell_y, 1, target] - 1
                tree_to_face[1, tree] = 2
            else # direction == 6, this is the +z face
                target = 1
                tree_to_tree[1, tree] = linear_indices[end - cell_y + 1, end,
                                                       target] - 1
                tree_to_face[1, tree] = 7 # first face dimensions are oppositely oriented, add 4
            end

            # Positive x-direction
            if cell_x < n_cells_x # Connect to tree at the same face
                tree_to_tree[2, tree] = linear_indices[cell_x + 1, cell_y,
                                                       direction] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 1 # This is the -x face
                target = 3
                tree_to_tree[2, tree] = linear_indices[1, cell_y, target] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 2 # This is the +x face
                target = 4
                tree_to_tree[2, tree] = linear_indices[1, cell_y, target] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 3 # This is the -y face
                target = 2
                tree_to_tree[2, tree] = linear_indices[1, cell_y, target] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 4 # This is the +y face
                target = 1
                tree_to_tree[2, tree] = linear_indices[1, cell_y, target] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 5 # This is the -z face
                target = 1
                tree_to_tree[2, tree] = linear_indices[end - cell_y + 1, 1,
                                                       target] - 1
                tree_to_face[2, tree] = 6 # first face dimensions are oppositely oriented, add 4
            else # direction == 6, this is the +z face
                target = 2
                tree_to_tree[2, tree] = linear_indices[cell_y, end, target] - 1
                tree_to_face[2, tree] = 3
            end

            # Negative y-direction
            if cell_y > 1 # Connect to tree at the same face
                tree_to_tree[3, tree] = linear_indices[cell_x, cell_y - 1,
                                                       direction] - 1
                tree_to_face[3, tree] = 3
            elseif direction == 1
                target = 5
                tree_to_tree[3, tree] = linear_indices[end, end - cell_x + 1,
                                                       target] - 1
                tree_to_face[3, tree] = 5 # first face dimensions are oppositely oriented, add 4
            elseif direction == 2
                target = 5
                tree_to_tree[3, tree] = linear_indices[1, cell_x, target] - 1
                tree_to_face[3, tree] = 0
            elseif direction == 3
                target = 5
                tree_to_tree[3, tree] = linear_indices[end - cell_x + 1, 1,
                                                       target] - 1
                tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4
            elseif direction == 4
                target = 5
                tree_to_tree[3, tree] = linear_indices[cell_x, end, target] - 1
                tree_to_face[3, tree] = 3
            elseif direction == 5
                target = 3
                tree_to_tree[3, tree] = linear_indices[end - cell_x + 1, 1,
                                                       target] - 1
                tree_to_face[3, tree] = 6 # first face dimensions are oppositely oriented, add 4
            else # direction == 6
                target = 3
                tree_to_tree[3, tree] = linear_indices[cell_x, end, target] - 1
                tree_to_face[3, tree] = 3
            end

            # Positive y-direction
            if cell_y < n_cells_y # Connect to tree at the same face
                tree_to_tree[4, tree] = linear_indices[cell_x, cell_y + 1,
                                                       direction] - 1
                tree_to_face[4, tree] = 2
            elseif direction == 1
                target = 6
                tree_to_tree[4, tree] = linear_indices[1, end - cell_x + 1,
                                                       target] - 1
                tree_to_face[4, tree] = 4 # first face dimensions are oppositely oriented, add 4
            elseif direction == 2
                target = 6
                tree_to_tree[4, tree] = linear_indices[end, cell_x, target] - 1
                tree_to_face[4, tree] = 1
            elseif direction == 3
                target = 6
                tree_to_tree[4, tree] = linear_indices[cell_x, 1, target] - 1
                tree_to_face[4, tree] = 2
            elseif direction == 4
                target = 6
                tree_to_tree[4, tree] = linear_indices[end - cell_x + 1, end,
                                                       target] - 1
                tree_to_face[4, tree] = 7 # first face dimensions are oppositely oriented, add 4
            elseif direction == 5
                target = 4
                tree_to_tree[4, tree] = linear_indices[cell_x, 1, target] - 1
                tree_to_face[4, tree] = 2
            else # direction == 6
                target = 4
                tree_to_tree[4, tree] = linear_indices[end - cell_x + 1, end,
                                                       target] - 1
                tree_to_face[4, tree] = 7 # first face dimensions are oppositely oriented, add 4
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

# Calculate physical coordinates of each node of a 2D cubed sphere mesh by mapping 
# the LGL quadrature nodes onto a Cartesian element grid, then using the equiangular 
# gnomonic mapping to place the nodes on the spherical surface
function calc_tree_node_coordinates_cubed_sphere_standard!(node_coordinates::AbstractArray{<:Any,
                                                                                           4},
                                                           nodes,
                                                           trees_per_face_dimension,
                                                           radius)
    n_cells_x = n_cells_y = trees_per_face_dimension

    linear_indices = LinearIndices((n_cells_x, n_cells_y, 6))

    # Get cell length in reference mesh
    dx = 2 / n_cells_x
    dy = 2 / n_cells_y

    for direction in 1:6
        for cell_y in 1:n_cells_y, cell_x in 1:n_cells_x
            tree = linear_indices[cell_x, cell_y, direction]

            x_offset = -1 + (cell_x - 1) * dx + dx / 2
            y_offset = -1 + (cell_y - 1) * dy + dy / 2
            z_offset = 0

            for j in eachindex(nodes), i in eachindex(nodes)
                # node_coordinates are the mapped reference node coordinates
                node_coordinates[:, i, j, tree] .= Trixi.cubed_sphere_mapping(x_offset +
                                                                              dx / 2 *
                                                                              nodes[i],
                                                                              y_offset +
                                                                              dy / 2 *
                                                                              nodes[j],
                                                                              z_offset,
                                                                              radius,
                                                                              0,
                                                                              direction)
            end
        end
    end
end

# Calculate physical coordinates of each node of a 2D cubed sphere mesh using the
# element-local mapping from Guba et al. (see https://doi.org/10.5194/gmd-7-2803-2014,
# Appendix A). We thank Oswald Knoth for bringing this paper to our attention.
function calc_tree_node_coordinates_cubed_sphere_local!(node_coordinates::AbstractArray{<:Any,
                                                                                        4},
                                                        nodes, trees_per_face_dimension,
                                                        radius)
    n_cells_x = n_cells_y = trees_per_face_dimension

    linear_indices = LinearIndices((n_cells_x, n_cells_y, 6))

    # Get cell length in reference mesh
    dx = 2 / n_cells_x
    dy = 2 / n_cells_y

    for direction in 1:6
        for cell_y in 1:n_cells_y, cell_x in 1:n_cells_x
            tree = linear_indices[cell_x, cell_y, direction]

            x_offset = -1 + (cell_x - 1) * dx + dx / 2
            y_offset = -1 + (cell_y - 1) * dy + dy / 2
            z_offset = 0

            # Vertices for bilinear mapping
            v1 = Trixi.cubed_sphere_mapping(x_offset - dx / 2,
                                            y_offset - dx / 2,
                                            z_offset, radius, 0, direction)

            v2 = Trixi.cubed_sphere_mapping(x_offset + dx / 2,
                                            y_offset - dx / 2,
                                            z_offset, radius, 0, direction)

            v3 = Trixi.cubed_sphere_mapping(x_offset + dx / 2,
                                            y_offset + dx / 2,
                                            z_offset, radius, 0, direction)

            v4 = Trixi.cubed_sphere_mapping(x_offset - dx / 2,
                                            y_offset + dx / 2,
                                            z_offset, radius, 0, direction)

            for j in eachindex(nodes), i in eachindex(nodes)
                # node_coordinates are the mapped reference node coordinates
                node_coordinates[:, i, j, tree] .= local_mapping(nodes[i], nodes[j],
                                                                 v1, v2, v3, v4, radius)
            end
        end
    end
end

# Local mapping from the reference quadrilateral to a spherical patch based on four vertex
# positions on the sphere, provided in Cartesian coordinates
@inline function local_mapping(xi1, xi2, v1, v2, v3, v4, radius)

    # Construct a bilinear mapping based on the four corner vertices
    x_bilinear = 0.25f0 * ((1 - xi1) * (1 - xi2) * v1 + (1 + xi1) * (1 - xi2) * v2 +
                  (1 + xi1) * (1 + xi2) * v3 + (1 - xi1) * (1 + xi2) * v4)

    # Project the mapped local coordinates onto the sphere
    return radius * x_bilinear / norm(x_bilinear)
end
end # @muladd
