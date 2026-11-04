# Utilities to complete `p4est` connectivities that only contain face connections.
#
# The connectivities of the meshes in this package are built by specifying the face
# neighbors of each tree (`tree_to_tree` and `tree_to_face`). However, `p4est` assumes
# that a connectivity is complete, i.e., that it also contains all corners (and edges in
# 3D) connecting trees that are not face neighbors. Otherwise, a single call to
# `p4est_balance` (`p8est_balance`) does not always produce a balanced forest, even for
# face balance (`P4EST_CONNECT_FACE`) as used by Trixi.jl. See
# https://github.com/cburstedde/p4est/issues/112 and
# https://github.com/trixi-framework/Trixi.jl/issues/1358.
#
# `p4est_connectivity_complete` computes the full connectivity from the vertices shared by
# the trees. Here, we derive the topological vertices from the face connections, so that
# the face connections (including their orientations) remain exactly as specified.

# Corners of each face in the `p4est` numbering, zero-based, cf. `p4est_face_corners`
# and `p8est_face_corners`
const P4EST_FACE_CORNERS = ((0, 2), (1, 3), (0, 1), (2, 3))
const P8EST_FACE_CORNERS = ((0, 2, 4, 6), (1, 3, 5, 7), (0, 1, 4, 5), (2, 3, 6, 7),
                            (0, 1, 2, 3), (4, 5, 6, 7))

face_corners_p4est(::Val{2}) = P4EST_FACE_CORNERS
face_corners_p4est(::Val{3}) = P8EST_FACE_CORNERS

# Number of the corner `corner` of face `face` as seen from the tree on the other side of
# this face, which is connected via its face `neighbor_face` with the given orientation
function face_neighbor_corner_p4est(::Val{2}, corner, face, neighbor_face, orientation)
    return Trixi.p4est_connectivity_face_neighbor_corner(corner, face, neighbor_face,
                                                         orientation)
end
function face_neighbor_corner_p4est(::Val{3}, corner, face, neighbor_face, orientation)
    return Trixi.p8est_connectivity_face_neighbor_corner(corner, face, neighbor_face,
                                                         orientation)
end

# Root of the set containing `i` in the union-find array `parent`, with path halving
function find_root!(parent, i)
    while parent[i] != i
        parent[i] = parent[parent[i]]
        i = parent[i]
    end
    return i
end

# Compute the topological vertices of a forest from the face connections of its trees.
# The corners of all trees are identified across all faces connecting trees. Since the
# trees around a vertex are connected via faces (also at physical boundaries), each
# resulting equivalence class of tree corners is a vertex of the forest.
# Return the zero-based `tree_to_vertex` array in `p4est` format and the number of
# vertices.
function vertices_from_faces(tree_to_tree, tree_to_face, ndims_::Val{NDIMS}) where {NDIMS}
    n_corners_per_tree = 2^NDIMS
    n_faces = 2 * NDIMS
    n_trees = size(tree_to_tree, 2)

    # Index of a tree corner in `parent`, with zero-based `tree` and `corner`
    corner_index(tree, corner) = tree * n_corners_per_tree + corner + 1

    parent = collect(1:(n_corners_per_tree * n_trees))
    for tree in 0:(n_trees - 1), face in 0:(n_faces - 1)
        neighbor = tree_to_tree[face + 1, tree + 1]
        # `tree_to_face` encodes the neighbor face and the orientation as
        # `neighbor_face + n_faces * orientation`
        orientation, neighbor_face = divrem(tree_to_face[face + 1, tree + 1], n_faces)

        # Physical boundary, the tree and face point to themselves
        if neighbor == tree && neighbor_face == face
            continue
        end

        for corner in face_corners_p4est(ndims_)[face + 1]
            neighbor_corner = face_neighbor_corner_p4est(ndims_, corner, face,
                                                         neighbor_face, orientation)
            root = find_root!(parent, corner_index(tree, corner))
            neighbor_root = find_root!(parent, corner_index(neighbor, neighbor_corner))
            parent[max(root, neighbor_root)] = min(root, neighbor_root)
        end
    end

    # Number the vertices consecutively in the order of their first appearance
    tree_to_vertex = Matrix{Trixi.p4est_topidx_t}(undef, n_corners_per_tree, n_trees)
    vertex_ids = Dict{Int, Trixi.p4est_topidx_t}()
    for tree in 0:(n_trees - 1), corner in 0:(n_corners_per_tree - 1)
        root = find_root!(parent, corner_index(tree, corner))
        tree_to_vertex[corner + 1, tree + 1] = get!(vertex_ids, root,
                                                    length(vertex_ids))
    end

    # `p4est` cannot represent trees touching themselves at different corners
    for tree in 1:n_trees
        if !allunique(view(tree_to_vertex, :, tree))
            error("tree $(tree - 1) touches itself, its connectivity cannot be completed")
        end
    end

    return tree_to_vertex, length(vertex_ids)
end

"""
    complete_connectivity(connectivity)

Return a complete copy of a `p4est` (`p8est`) `connectivity` that only contains face
connections, i.e., a connectivity that additionally contains all corners (and edges in 3D)
connecting trees. The face connections remain unchanged. The original `connectivity` is
destroyed.

This is required for `p4est_balance` (`p8est_balance`) to produce a balanced forest, see
https://github.com/cburstedde/p4est/issues/112.
"""
function complete_connectivity(connectivity::Ptr{Trixi.p4est_connectivity_t})
    connectivity_pw = Trixi.PointerWrapper(connectivity)
    n_trees = Int(connectivity_pw.num_trees[])
    tree_to_tree = copy(unsafe_wrap(Array, connectivity_pw.tree_to_tree, (4, n_trees)))
    tree_to_face = copy(unsafe_wrap(Array, connectivity_pw.tree_to_face, (4, n_trees)))

    tree_to_vertex, n_vertices = vertices_from_faces(tree_to_tree, tree_to_face, Val(2))
    # `p4est_connectivity_complete` only uses the topology given by `tree_to_vertex`.
    # The vertex coordinates are only used by `p4est` to write VTK files, which Trixi.jl
    # does not use, so they can be zero.
    vertices = zeros(Float64, 3, n_vertices)

    # The corner connectivity is computed by `p4est_connectivity_complete` below
    n_corners = 0
    tree_to_corner = C_NULL
    ctt_offset = zeros(Trixi.p4est_topidx_t, 1)
    corner_to_tree = C_NULL
    corner_to_corner = C_NULL

    connectivity_complete = Trixi.p4est_connectivity_new_copy(n_vertices, n_trees,
                                                              n_corners, vertices,
                                                              tree_to_vertex,
                                                              tree_to_tree, tree_to_face,
                                                              tree_to_corner, ctt_offset,
                                                              corner_to_tree,
                                                              corner_to_corner)
    Trixi.p4est_connectivity_complete(connectivity_complete)
    Trixi.p4est_connectivity_destroy(connectivity)

    @assert Trixi.p4est_connectivity_is_valid(connectivity_complete) == 1
    check_face_connections(connectivity_complete, tree_to_tree, tree_to_face)

    return connectivity_complete
end

function complete_connectivity(connectivity::Ptr{Trixi.p8est_connectivity_t})
    connectivity_pw = Trixi.PointerWrapper(connectivity)
    n_trees = Int(connectivity_pw.num_trees[])
    tree_to_tree = copy(unsafe_wrap(Array, connectivity_pw.tree_to_tree, (6, n_trees)))
    tree_to_face = copy(unsafe_wrap(Array, connectivity_pw.tree_to_face, (6, n_trees)))

    tree_to_vertex, n_vertices = vertices_from_faces(tree_to_tree, tree_to_face, Val(3))
    # `p8est_connectivity_complete` only uses the topology given by `tree_to_vertex`.
    # The vertex coordinates are only used by `p4est` to write VTK files, which Trixi.jl
    # does not use, so they can be zero.
    vertices = zeros(Float64, 3, n_vertices)

    # The edge and corner connectivity is computed by `p8est_connectivity_complete` below
    n_edges = 0
    tree_to_edge = C_NULL
    ett_offset = zeros(Trixi.p4est_topidx_t, 1)
    edge_to_tree = C_NULL
    edge_to_edge = C_NULL

    n_corners = 0
    tree_to_corner = C_NULL
    ctt_offset = zeros(Trixi.p4est_topidx_t, 1)
    corner_to_tree = C_NULL
    corner_to_corner = C_NULL

    connectivity_complete = Trixi.p8est_connectivity_new_copy(n_vertices, n_trees,
                                                              n_edges, n_corners,
                                                              vertices, tree_to_vertex,
                                                              tree_to_tree, tree_to_face,
                                                              tree_to_edge, ett_offset,
                                                              edge_to_tree, edge_to_edge,
                                                              tree_to_corner, ctt_offset,
                                                              corner_to_tree,
                                                              corner_to_corner)
    Trixi.p8est_connectivity_complete(connectivity_complete)
    Trixi.p8est_connectivity_destroy(connectivity)

    @assert Trixi.p8est_connectivity_is_valid(connectivity_complete) == 1
    check_face_connections(connectivity_complete, tree_to_tree, tree_to_face)

    return connectivity_complete
end

# `p4est_connectivity_complete` recomputes the face connections from the vertices. Make
# sure that they agree with the ones specified originally, including the orientations.
function check_face_connections(connectivity, tree_to_tree, tree_to_face)
    connectivity_pw = Trixi.PointerWrapper(connectivity)
    tree_to_tree_new = unsafe_wrap(Array, connectivity_pw.tree_to_tree,
                                   size(tree_to_tree))
    tree_to_face_new = unsafe_wrap(Array, connectivity_pw.tree_to_face,
                                   size(tree_to_face))
    if tree_to_tree_new != tree_to_tree || tree_to_face_new != tree_to_face
        error("completing the connectivity changed the face connections")
    end

    return nothing
end
