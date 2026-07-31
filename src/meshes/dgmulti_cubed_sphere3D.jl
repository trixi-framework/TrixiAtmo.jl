function DGMultiMeshCubedShell3D(dg::DGMulti{3, <:Hex}, inner_radius, outer_radius;
                                       horizontal_initial_refinement = 3,
                                       vertical_layers = 5)
    NDIMS_AMBIENT = 3

    radius = 1.0

    vertex_coordinates = calc_node_coordinates_cube_vertices(radius)

    Vxyz_quad = ntuple(n -> vertex_coordinates[n, :], NDIMS_AMBIENT)

    EToV_quad = zeros(Int, 6, 4)

    for i in 1:size(EToV_quad, 1)
        EToV_quad[i, :] = cubed_quadrilaterals_vertices_idx_map[i]
    end

    for j in 1:horizontal_initial_refinement
        EToV_old = EToV_quad
        Vxyz_old = ntuple(n -> copy(Vxyz_quad[n]), NDIMS_AMBIENT)
        old_to_new = Dict{Int, Int}()
        edge_to_new = Dict{Tuple{Int, Int}, Int}()
        EToV_quad = zeros(Int, (size(EToV_old, 1) * 4, 4))

        Vxyz_quad = ntuple(n -> Vector{eltype(Vxyz_quad[1])}(), NDIMS_AMBIENT)

        for i in 1:size(EToV_old, 1)
            idx_old = EToV_old[i, :]

            for k in idx_old
                if !haskey(old_to_new, k)
                    old_to_new[k] = length(Vxyz_quad[1]) + 1
                    for n in 1:NDIMS_AMBIENT
                        push!(Vxyz_quad[n], Vxyz_old[n][k])
                    end
                end
            end

            k_counter = 1
            for k in idx_old
                l_counter =  1
                for l in idx_old
                    if k < l && l_counter + k_counter != 5
                        edge = (k, l)
                        if !haskey(edge_to_new, edge)
                            edge_to_new[edge] = length(Vxyz_quad[1]) + 1
                            vk = ntuple(n -> Vxyz_old[n][k], NDIMS_AMBIENT)
                            vl = ntuple(n -> Vxyz_old[n][l], NDIMS_AMBIENT)
                            midpoint = 0.5 .* (vk .+ vl)
                            midpoint = midpoint .* (radius / norm(midpoint)) # Normalize to outer radius
                            for n in 1:NDIMS_AMBIENT
                                push!(Vxyz_quad[n], midpoint[n])
                            end
                        end
                    end
                    l_counter += 1
                end
                k_counter += 1
            end

            midpoint = ntuple(n -> 0.25 * sum(Vxyz_old[n][idx_old]), NDIMS_AMBIENT)
            midpoint = midpoint .* (radius / norm(midpoint))                        # Normalize to outer radius
            for n in 1:NDIMS_AMBIENT
                push!(Vxyz_quad[n], midpoint[n])
            end

            id1 = old_to_new[idx_old[1]]
            id2 = edge_to_new[Tuple(sort([idx_old[1], idx_old[2]]))]
            id3 = old_to_new[idx_old[2]]
            id4 = edge_to_new[Tuple(sort([idx_old[1], idx_old[3]]))]
            id5 = length(Vxyz_quad[1])
            id6 = edge_to_new[Tuple(sort([idx_old[2], idx_old[4]]))]
            id7 = old_to_new[idx_old[3]]
            id8 = edge_to_new[Tuple(sort([idx_old[3], idx_old[4]]))]
            id9 = old_to_new[idx_old[4]]

            ids = [id1, id2, id3, id4, id5, id6, id7, id8, id9]

            # Fill EToV for the 4 new triangles
            for (sub_i, vertex_map) in enumerate(cubed_quad_vertices_idx_map)
                EToV_quad[(i - 1) * 4 + sub_i, :] = getindex.(Ref(ids), vertex_map)
            end

        end
    end

    # Create the prism elements by extruding the 2D mesh
    num_quads = size(EToV_quad, 1)
    num_points_per_layer = size(Vxyz_quad[1], 1)
    EToV = zeros(Int, num_quads * vertical_layers, 8)
    Vxyz = ntuple(n -> Vxyz_quad[n] .* inner_radius, NDIMS_AMBIENT)
    radius_list = range(inner_radius, outer_radius, length=vertical_layers + 1)

    for i in 1:vertical_layers
        EToV[(i - 1) * num_quads + 1:i * num_quads, 1:4] = EToV_quad .+ (i - 1) * num_points_per_layer
        EToV[(i - 1) * num_quads + 1:i * num_quads, 5:8] = EToV_quad .+ i * num_points_per_layer
        Vxyz_layer = ntuple(n -> Vxyz_quad[n] .* radius_list[i + 1], NDIMS_AMBIENT)
        Vxyz = ntuple(n -> vcat(Vxyz[n], Vxyz_layer[n]), NDIMS_AMBIENT)
    end
    md = StartUpDG.MeshData(Vxyz, EToV, dg.basis)
    project_onto_sphere!(md, dg)
    on_bottom_boundary((x, y, z), tol = 1e-13) = x^2 + y^2 + z^2 - inner_radius^2 < tol
    on_top_boundary((x, y, z), tol = 1e-13) = x^2 + y^2 + z^2 - outer_radius^2 < tol
    boundary_faces = StartUpDG.tag_boundary_faces(md, Dict(:bottom => on_bottom_boundary,
                                                       :top    => on_top_boundary))
    return DGMultiMesh(dg, Trixi.GeometricTermsType(Trixi.Curved(), dg), md, boundary_faces)
end

function project_onto_sphere!(md::MeshData, dg::DGMulti{NDIMS, <:Hex}) where {NDIMS}
    rd = dg.basis
    (; xyz, xyzq, xyzf) = md
    for e in 1:size(md.xyz[1], 2)
        VX, VY, VZ = map(coords -> transpose(coords[:, e]) / rd.V1', md.xyz)
        vertices = ntuple(n -> [VX[n], VY[n], VZ[n]], 8)
        inner_radius = norm(vertices[1])
        outer_radius = norm(vertices[5])
        for j in 1:size(rd.rst[1], 1)
            t = rd.rst[3][j]
            radius = (1 - t) / 2 * inner_radius + (1 + t) / 2 * outer_radius
            x_node = ntuple(n -> xyz[n][j, e], NDIMS)
            x_node = radius / norm(x_node) .* x_node

            for n in 1:3
                xyz[n][j, e] = x_node[n]
            end
        end

        for j in 1:size(rd.rstq[1], 1)
            r, s, t = rd.rstq[1][j], rd.rstq[2][j], rd.rstq[3][j]
            radius = (1 - t) / 2 * inner_radius + (1 + t) / 2 * outer_radius
            x_node = ntuple(n -> xyzq[n][j, e], NDIMS)
            x_node = radius / norm(x_node) .* x_node

            for n in 1:3
                xyzq[n][j, e] = x_node[n]
            end
        end

        for j in 1:size(rd.rstf[1], 1)
            r, s, t = rd.rstf[1][j], rd.rstf[2][j], rd.rstf[3][j]
            radius = (1 - t) / 2 * inner_radius + (1 + t) / 2 * outer_radius
            x_node = ntuple(n -> xyzf[n][j, e], NDIMS)
            x_node = radius / norm(x_node) .* x_node

            for n in 1:3
                xyzf[n][j, e] = x_node[n]
            end
        end
    end
    md = setproperties(md, xyz = xyz, xyzq = xyzq, xyzf = xyzf)
    return md
end

function calc_node_coordinates_cube_vertices(radius; RealT = Float64)
    vertices = Array{RealT, 2}(undef, 3, 12)

    vertices[:, 1] = [ 1/sqrt(3),  1/sqrt(3),  1/sqrt(3)]
    vertices[:, 2] = [-1/sqrt(3),  1/sqrt(3),  1/sqrt(3)]
    vertices[:, 3] = [ 1/sqrt(3), -1/sqrt(3),  1/sqrt(3)]
    vertices[:, 4] = [-1/sqrt(3), -1/sqrt(3),  1/sqrt(3)]
    vertices[:, 5] = [ 1/sqrt(3),  1/sqrt(3), -1/sqrt(3)]
    vertices[:, 6] = [-1/sqrt(3),  1/sqrt(3), -1/sqrt(3)]
    vertices[:, 7] = [ 1/sqrt(3), -1/sqrt(3), -1/sqrt(3)]
    vertices[:, 8] = [-1/sqrt(3), -1/sqrt(3), -1/sqrt(3)]

    return vertices * radius
end

# Index map for the vertices of each triangle on the triangular faces of the icosahedron (see Fig 4)
const cubed_quad_vertices_idx_map = ([1, 2, 4, 5], # Quad 1
                                     [2, 3, 5, 6], # Tri 2
                                     [4, 5, 7, 8], # Tri 3
                                     [5, 6, 8, 9]) # Tri 4

const cubed_quadrilaterals_vertices_idx_map = ([1, 2, 3, 4],
                                               [1, 5, 2, 6],
                                               [1, 3, 5, 7],
                                               [3, 4, 7, 8],
                                               [2, 4, 6, 8],
                                               [5, 6, 7, 8])
                                               