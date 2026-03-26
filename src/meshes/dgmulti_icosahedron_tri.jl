"""
    DGMultiMeshTriIcosahedron2D(trees_per_face_dimension, radius;
                               polydeg, RealT=Float64,
                               initial_refinement_level=0, unsaved_changes=true,
                               p4est_partition_allow_for_coarsening=true)

Build a triangle-based icosahedral mesh as a 2D `DGMultiMesh` with 20 unrefined triangular faces. The mesh
is then refined uniformly to the specified `initial_refinement_level`, yielding a total of `20 * 4^initial_refinement_level` triangles.

The mesh will have no boundaries.

# Arguments
- `dg::DGMulti{NDIMS, <:Tri}`: the DGMulti discretization to use for the mesh.
- `radius::Integer`: the radius of the sphere.
- `initial_refinement_level::Integer`: refine the mesh uniformly to this level before the
  simulation starts.
"""
function DGMultiMeshTriIcosahedron2D(dg::DGMulti{2, <:Tri}, radius;
                                     initial_refinement_level = 3,
                                     is_on_boundary = nothing)
    NDIMS_AMBIENT = 3

    vertex_coordinates = calc_node_coordinates_icosahedron_vertices(radius)

    Vxyz = ntuple(n -> vertex_coordinates[n, :], NDIMS_AMBIENT)

    EToV = zeros(Int, 20, 3)

    for i in 1:size(EToV, 1)
        EToV[i, :] = icosahedron_triangle_vertices_idx_map[i]
    end

    for j in 1:initial_refinement_level
        EToV_old = EToV
        Vxyz_old = ntuple(n -> copy(Vxyz[n]), NDIMS_AMBIENT)
        # We need to keep track of the old vertex indices and edge indices to avoid creating duplicate vertices when refining the mesh.
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
                        midpoint = midpoint .* (radius / norm(midpoint)) # Normalize to outer radius
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
    md = project_onto_sphere!(md, dg, radius)
    boundary_faces = StartUpDG.tag_boundary_faces(md, is_on_boundary)
    return DGMultiMesh(dg, Trixi.GeometricTermsType(Trixi.Curved(), dg), md, boundary_faces)
end

function project_onto_sphere!(md::MeshData, dg::DGMulti{NDIMS}, radius) where {NDIMS}
    x, y, z = md.xyz
    norms = sqrt.(x.^2 + y.^2 + z.^2)
    xyz = ntuple(n -> md.xyz[n] ./ norms .* radius, 3)
    xq, yq, zq = md.xyzq
    norms_q = sqrt.(xq.^2 + yq.^2 + zq.^2)
    xyzq = ntuple(n -> md.xyzq[n] ./ norms_q .* radius, 3)
    xf, yf, zf = md.xyzf
    norms_f = sqrt.(xf.^2 + yf.^2 + zf.^2)
    xyzf = ntuple(n -> md.xyzf[n] ./ norms_f .* radius, 3)
    md = setproperties(md, xyz = xyz, xyzq = xyzq, xyzf = xyzf)
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

function StartUpDG.geometric_factors(x, y, z, Dr, Ds; Filters = (I, I, I))
    xr, xs = Dr * x, Ds * x
    yr, ys = Dr * y, Ds * y
    zr, zs = Dr * z, Ds * z

    a1 = Dr * x, Dr * y, Dr * z
    a2 = Ds * x, Dr * y, Ds * z

    a3 = @. a1[2] * a2[3] - a1[3] * a2[2],
            a1[3] * a2[1] - a1[1] * a2[3],
            a1[1] * a2[2] - a1[2] * a2[1]

    norm_a3 = sqrt.(a3[1] .^ 2 + a3[2] .^ 2 + a3[3] .^ 2)

    a3 = @. a3[1] / norm_a3, a3[2] / norm_a3, a3[3] / norm_a3

    rxJ, ryJ, rzJ = @. a2[2] * a3[3] - a2[3] * a3[2],
                       a2[3] * a3[1] - a2[1] * a3[3],
                       a2[1] * a3[2] - a2[2] * a3[1]
    sxJ, syJ, szJ = @. a3[2] * a1[3] - a3[3] * a1[2],
                       a3[3] * a1[1] - a3[1] * a1[3],
                       a3[1] * a1[2] - a3[2] * a1[1]

    J = @. sqrt((xr^2 + yr^2 + zr^2) * (xs^2 + ys^2 + zs^2) -
                (xr * xs + yr * ys + zr * zs)^2)
    return rxJ, sxJ, ryJ, syJ, rzJ, szJ, J
end

# physical outward normals are computed via Nanson's formula: G * nhatJ, where 
# G = matrix of J-scaled geometric terms. Here, Vf is a face interpolation matrix 
# which maps interpolation nodes to face nodes. 
function StartUpDG.compute_normals(geo::SMatrix{Dim, DimAmbient}, Vf,
                                   nrstJ...) where {Dim, DimAmbient}
    nxyzJ = ntuple(x -> zeros(size(Vf, 1), size(first(geo), 2)), DimAmbient)
    for i in 1:Dim, j in 1:DimAmbient
        nxyzJ[i] .+= (Vf * geo[i, j]) .* nrstJ[j]
    end
    Jf = sqrt.(sum(map(x -> x .^ 2, nxyzJ)))
    return nxyzJ..., Jf
end

function StartUpDG.MeshData(VX, VY, VZ, EToV, rd::RefElemData{2};
                            is_periodic = (false, false, false))
    (; fv) = rd
    FToF = StartUpDG.connect_mesh(EToV, fv)
    Nfaces, K = size(FToF)

    #Construct global coordinates
    (; V1) = rd
    x, y, z = (x -> V1 * x[transpose(EToV)]).((VX, VY, VZ))

    #Compute connectivity maps: uP = exterior value used in DG numerical fluxes
    (; r, s, Vf) = rd
    xf, yf, zf = (x -> Vf * x).((x, y, z))
    mapM, mapP, mapB = StartUpDG.build_node_maps(FToF, (xf, yf, zf))
    mapM = reshape(mapM, :, K)
    mapP = reshape(mapP, :, K)

    #Compute geometric factors and surface normals
    (; Dr, Ds) = rd
    rxJ, sxJ, ryJ, syJ, rzJ, szJ, J = StartUpDG.geometric_factors(x, y, z, Dr, Ds)
    rstxyzJ = SMatrix{2, 3}(rxJ, ryJ, rzJ, sxJ, syJ, szJ)
    nrstJ = (rd.nrstJ..., zeros(size(rd.nrstJ[1])))
    (; Vq, wq) = rd
    xq, yq, zq = (x -> Vq * x).((x, y, z))
    wJq = diagm(wq) * (Vq * J)

    nxJ, nyJ, nzJ, Jf = StartUpDG.compute_normals(rstxyzJ, rd.Vf, nrstJ...)

    periodicity = (false, false, false)
    md = MeshData(StartUpDG.VertexMappedMesh(rd.element_type, tuple(VX, VY, VZ), EToV),
                  FToF, tuple(x, y, z), tuple(xf, yf, zf), tuple(xq, yq, zq), wJq,
                  mapM, mapP, mapB,
                  rstxyzJ, J, tuple(nxJ, nyJ, nzJ), Jf,
                  periodicity)

    if any(is_periodic)
        # loosen the tolerance if N >> 1
        tol = length(rd.r) * 100 * eps()
        md = make_periodic(md, is_periodic; tol)
    end

    return md
end
