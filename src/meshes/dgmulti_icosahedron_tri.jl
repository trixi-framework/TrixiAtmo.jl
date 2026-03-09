function DGMultiMeshTriIcosahedron2D(dg::DGMulti{NDIMS, <:Tri}, radius;
                                     initial_refinement = 3,
                                     is_on_boundary = nothing) where {NDIMS}

    NDIMS_AMBIENT = 3

    vertex_coordinates = calc_node_coordinates_icosahedron_vertices(radius)

    Vxyz = ntuple(n -> vertex_coordinates[n, :], NDIMS_AMBIENT)

    EToV = zeros(Int, 20, 3)
    
    for i in 1:size(EToV, 1)
        EToV[i, :] = icosahedron_triangle_vertices_idx_map[i]
    end

    for j in 1:initial_refinement
        EToV_old = EToV
        Vxyz_old = ntuple(n -> copy(Vxyz[n]), NDIMS_AMBIENT)
        old_to_new = Dict{Int, Int}()
        edge_to_new = Dict{Tuple{Int, Int}, Int}()
        EToV = zeros(Int, (size(EToV_old, 1) * 4, 3))
        
        Vxyz = ntuple(n -> Vector{eltype(Vxyz[1])}(), NDIMS_AMBIENT)

        for i in 1:size(EToV_old, 1)
            idx_old = EToV_old[i, :]

            for k in idx_old
                if !haskey(old_to_new, k)
                    old_to_new[k] = length(Vxyz[1]) + 1
                    for n in 1:NDIMS_AMBIENT
                        push!(Vxyz[n], Vxyz_old[n][k])
                    end
                end
            end

            for k in idx_old, l in idx_old
                if k < l
                    edge = (k, l)
                    if !haskey(edge_to_new, edge)
                        edge_to_new[edge] = length(Vxyz[1]) + 1
                        vk = ntuple(n -> Vxyz_old[n][k], NDIMS_AMBIENT)
                        vl = ntuple(n -> Vxyz_old[n][l], NDIMS_AMBIENT)
                        midpoint = 0.5 .* (vk .+ vl)
                        midpoint = midpoint .* (EARTH_RADIUS / norm(midpoint)) # Normalize to outer radius
                        for n in 1:NDIMS_AMBIENT
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

    md = StartUpDG.MeshData(Vxyz, EToV, dg.basis)
    md = spherify_meshdata!(md, dg, EToV, Vxyz, radius)
    boundary_faces = StartUpDG.tag_boundary_faces(md, is_on_boundary)
    return DGMultiMesh(dg, Trixi.GeometricTermsType(Trixi.Curved(), dg), md, boundary_faces)
end

function spherify_meshdata!(md::MeshData, dg::DGMulti{NDIMS}, EToV, Vxyz, radius) where {NDIMS}
    rd = dg.basis
    (; xyz, xyzq, xyzf) = md
    for e in 1:size(EToV, 1)
        v1, v2, v3 = ntuple(n -> SVector(ntuple(d -> Vxyz[d][EToV[e, n]], 3)), 3)
        for j in 1:size(rd.rst[1], 1)
            r, s = rd.rst[1][j], rd.rst[2][j]
            # Bilinear mapping from reference square to physical space
            x_node = local_mapping(r, s, v1, v2, v3, radius)

            for n in 1:3
                xyz[n][j, e] = x_node[n]
            end
        end

        for j in 1:size(rd.rstq[1], 1)
            r, s = rd.rstq[1][j], rd.rstq[2][j]
            # Bilinear mapping from reference square to physical space
            x_node = local_mapping(r, s, v1, v2, v3, radius)

            for n in 1:3
                xyzq[n][j, e] = x_node[n]
            end
        end

        for j in 1:size(rd.rstf[1], 1)
            r, s = rd.rstf[1][j], rd.rstf[2][j]
            # Bilinear mapping from reference square to physical space
            x_node = local_mapping(r, s, v1, v2, v3, radius)

            for n in 1:3
                xyzf[n][j, e] = x_node[n]
            end
        end
    end
    md = @set md.xyz = xyz
    md = @set md.xyzq = xyzq
    md = @set md.xyzf = xyzf
    return md
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
#            /⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺⎺\
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

# Local mapping from the reference triangle to a spherical patch based on three vertex
# positions on the sphere, provided in Cartesian coordinates
@inline function local_mapping(r, s, v1, v2, v3, radius)

    # Construct a bilinear mapping based on the four corner vertices
    xe = 0.5f0 * (-(r + s) * v1 +
                   (1 + r) * v2 +
                   (1 + s) * v3)

    # Project the mapped local coordinates onto the sphere
    return radius * xe / norm(xe)
end