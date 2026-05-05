using Trixi: connectivity_cubed_sphere, new_p4est

abstract type AdaptVerticalGrid end

struct GalChen <: AdaptVerticalGrid end

struct Sleve{RealT} <: AdaptVerticalGrid
    etaH::RealT
    s::RealT
end

@inline function (adapt_galchen::GalChen)(z_reference, z_topography, H)
    z = z_reference + (H - z_reference) * z_topography / H
    return z
end

@inline function (adapt_sleve::Sleve)(z_reference, z_topography, H)
    (; s, etaH) = Sleve
    eta = z_reference / H
    if eta <= etaH
        z = eta * H +
            z_topography * sinh((etaH - eta) / s / etaH) / sinh(one(z_reference) / s)
    else
        z = eta * H
    end
    return z
end

"""
    P4estMeshCubedSphereTopography(trees_per_face_dimension, layers, inner_radius, thickness;
                         polydeg, RealT=Float64,
                         initial_refinement_level=0, unsaved_changes=true,
                         p4est_partition_allow_for_coarsening=true)

Build a "Cubed Sphere" mesh as `P4estMesh` with
`6 * trees_per_face_dimension^2 * layers` trees.

The mesh will have two boundaries, `:inside` and `:outside`.

# Arguments
- `trees_per_face_dimension::Integer`: the number of trees in the first two local dimensions of
                                       each face.
- `layers::Integer`: the number of trees in the third local dimension of each face, i.e., the number
                     of layers of the sphere.
- `inner_radius::RealT`: the inner radius of the sphere.
- `thickness::RealT`: the thickness of the spherical shell. The outer radius will be `inner_radius + thickness`.
- `polydeg::Integer`: polynomial degree used to store the geometry of the mesh.
                      The mapping will be approximated by an interpolation polynomial
                      of the specified degree for each tree.
- `RealT::Type`: the type that should be used for coordinates.
- `initial_refinement_level::Integer`: refine the mesh uniformly to this level before the simulation starts.
- `unsaved_changes::Bool`: if set to `true`, the mesh will be saved to a mesh file.
- `p4est_partition_allow_for_coarsening::Bool`: Must be `true` when using AMR to make mesh adaptivity
                                                independent of domain partitioning. Should be `false` for static meshes
                                                to permit more fine-grained partitioning.
"""
function P4estMeshCubedSphereTopography(trees_per_face_dimension, layers, inner_radius,
                                        thickness;
                                        polydeg, RealT = Float64,
                                        initial_refinement_level = 0,
                                        unsaved_changes = true,
                                        p4est_partition_allow_for_coarsening = true,
                                        initial_topography, adapt_vertical_grid = GalChen())
    connectivity = connectivity_cubed_sphere(trees_per_face_dimension, layers)

    n_trees = 6 * trees_per_face_dimension^2 * layers

    basis = LobattoLegendreBasis(RealT, polydeg)
    nodes = basis.nodes

    tree_node_coordinates = Array{RealT, 5}(undef, 3,
                                            ntuple(_ -> length(nodes), 3)...,
                                            n_trees)
    calc_tree_node_coordinates!(tree_node_coordinates, nodes,
                                trees_per_face_dimension, layers,
                                inner_radius, thickness, initial_topography,
                                adapt_vertical_grid)

    p4est = new_p4est(connectivity, initial_refinement_level)

    boundary_names = fill(Symbol("---"), 2 * 3, n_trees)
    boundary_names[5, :] .= Symbol("inside")
    boundary_names[6, :] .= Symbol("outside")

    return P4estMesh{3}(p4est, tree_node_coordinates, nodes,
                        boundary_names, "", unsaved_changes,
                        p4est_partition_allow_for_coarsening)
end

function cubed_sphere_mapping_topography(xi, eta, zeta, inner_radius, thickness, direction,
                                         initial_topography, adapt_vertical_grid)
    alpha = xi * pi / 4
    beta = eta * pi / 4

    x = tan(alpha)
    y = tan(beta)

    cube_coordinates = (SVector(-1, -x, y),
                        SVector(1, x, y),
                        SVector(x, -1, y),
                        SVector(-x, 1, y),
                        SVector(-x, y, -1),
                        SVector(x, y, 1))

    r = sqrt(1 + x^2 + y^2)

    unit_pt = cube_coordinates[direction] / r
    z_topography = initial_topography(unit_pt[1], unit_pt[2], unit_pt[3])

    z_reference = thickness * 0.5 * (zeta + 1)
    z = adapt_vertical_grid(z_reference, z_topography, thickness)
    return (inner_radius + z) / r * cube_coordinates[direction]
end

function calc_tree_node_coordinates!(node_coordinates::AbstractArray{<:Any, 5},
                                     nodes, trees_per_face_dimension, layers,
                                     inner_radius, thickness, initial_topography,
                                     adapt_vertical_grid)
    n_cells_x = n_cells_y = trees_per_face_dimension
    n_cells_z = layers

    linear_indices = LinearIndices((n_cells_x, n_cells_y, n_cells_z, 6))

    dx = 2 / n_cells_x
    dy = 2 / n_cells_y
    dz = 2 / n_cells_z

    for direction in 1:6
        for cell_z in 1:n_cells_z, cell_y in 1:n_cells_y, cell_x in 1:n_cells_x
            tree = linear_indices[cell_x, cell_y, cell_z, direction]

            x_offset = -1 + (cell_x - 1) * dx + dx / 2
            y_offset = -1 + (cell_y - 1) * dy + dy / 2
            z_offset = -1 + (cell_z - 1) * dz + dz / 2

            for k in eachindex(nodes), j in eachindex(nodes), i in eachindex(nodes)
                node_coordinates[:, i, j, k, tree] .= cubed_sphere_mapping_topography(x_offset +
                                                                                      dx /
                                                                                      2 *
                                                                                      nodes[i],
                                                                                      y_offset +
                                                                                      dy /
                                                                                      2 *
                                                                                      nodes[j],
                                                                                      z_offset +
                                                                                      dz /
                                                                                      2 *
                                                                                      nodes[k],
                                                                                      inner_radius,
                                                                                      thickness,
                                                                                      direction,
                                                                                      initial_topography,
                                                                                      adapt_vertical_grid)
            end
        end
    end

    return nothing
end
