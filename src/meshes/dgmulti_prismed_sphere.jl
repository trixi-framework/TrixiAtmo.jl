# Construct a DGMultiMesh consisting of prismatic elements based on a refined icosahedron
# As of now this 3D mesh only serves to simulate 2D covariant equations on the sphere, with
# the spherical shell thus only consisting of a single layer of prismatic elements extruded
# from the triangular faces of the icosahedron.
function DGMultiMeshPrismIcosahedron(dg::DGMulti{NDIMS};
    inner_radius = 1.0 * EARTH_RADIUS,
    outer_radius = 0.99 * EARTH_RADIUS,
    initial_refinement = 3,
    is_on_boundary = nothing) where {NDIMS}

    vertex_coordinates = calc_node_coordinates_icosahedron_vertices(outer_radius)

    Vxyz = ntuple(n -> vertex_coordinates[n, :], NDIMS)

    EToV = zeros(Int, 20, 3)
    
    for i in 1:size(EToV, 1)
        EToV[i, :] = icosahedron_triangle_vertices_idx_map[i]
    end

    for j in 1:initial_refinement
        EToV_old = EToV
        Vxyz_old = ntuple(n -> copy(Vxyz[n]), NDIMS)
        old_to_new = Dict{Int, Int}()
        edge_to_new = Dict{Tuple{Int, Int}, Int}()
        EToV = zeros(Int, (size(EToV_old, 1) * 4, 3))
        
        Vxyz = ntuple(n -> Vector{eltype(Vxyz[1])}(), NDIMS)

        for i in 1:size(EToV_old, 1)
            idx_old = EToV_old[i, :]

            for k in idx_old
                if !haskey(old_to_new, k)
                    old_to_new[k] = length(Vxyz[1]) + 1
                    for n in 1:NDIMS
                        push!(Vxyz[n], Vxyz_old[n][k])
                    end
                end
            end

            for k in idx_old, l in idx_old
                if k < l
                    edge = (k, l)
                    if !haskey(edge_to_new, edge)
                        edge_to_new[edge] = length(Vxyz[1]) + 1
                        vk = ntuple(n -> Vxyz_old[n][k], NDIMS)
                        vl = ntuple(n -> Vxyz_old[n][l], NDIMS)
                        midpoint = 0.5 .* (vk .+ vl)
                        midpoint = midpoint .* (outer_radius / norm(midpoint)) # Normalize to outer radius
                        for n in 1:NDIMS
                            push!(Vxyz[n], midpoint[n])
                        end
                    end
                end
            end
            id1 = old_to_new[idx_old[1]]
            id2 = edge_to_new[Tuple(sort([idx_old[1], idx_old[2]]))]
            id3 = old_to_new[idx_old[2]]
            id4 = edge_to_new[Tuple(sort([idx_old[2], idx_old[3]]))]
            id5 = old_to_new[idx_old[3]]
            id6 = edge_to_new[Tuple(sort([idx_old[3], idx_old[1]]))]

            ids = [id1, id2, id3, id4, id5, id6]

            # Fill EToV for the 4 new triangles
            for (sub_i, vertex_map) in enumerate(icosahedron_tri_vertices_idx_map)
                EToV[(i - 1) * 4 + sub_i, :] = getindex.(Ref(ids), vertex_map)
            end

        end
    end

    num_vertices = size(Vxyz[1], 1)
    vertices_inner = ntuple(n -> Vxyz[n] .* (inner_radius / outer_radius), NDIMS)

    Vxyz = ntuple(n -> vcat(Vxyz[n], vertices_inner[n]), NDIMS)
    EToV = hcat(EToV .+ num_vertices, EToV)

    md = MeshData(Vxyz, EToV, dg.basis)
    boundary_faces = StartUpDG.tag_boundary_faces(md, is_on_boundary)
    return DGMultiMesh(dg, Trixi.GeometricTermsType(Trixi.Curved(), dg), md, boundary_faces)
end

# We use a local numbering to obtain the triangle vertices of each triangular face
#
# Fig: Local quad vertex numbering for a triangular face (corner vertices of the triangular face in parenthesis)
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
#            /⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺\
#           / \                 /  \
#          /   \               /    \
#         /     \             /      \
#        /       \           /        \
#       /         \         /          \ 
#      /           \       /            \
#     /             \     /              \
#    /               \   /                \
#  1/_________________\_/_________________3\
#   (1)                2                    (2)

# Index map for the vertices of each triangle on the triangular faces of the icosahedron (see Fig 4)
const icosahedron_tri_vertices_idx_map = ([1, 2, 6], # Tri 1
                                          [2, 3, 4], # Tri 2
                                          [4, 5, 6], # Tri 3
                                          [2, 4, 6]) # Tri 4