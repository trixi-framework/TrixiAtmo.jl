using Trixi: new_p4est, p4est_topidx_t, p8est_connectivity_new_copy,
             p8est_connectivity_is_valid

abstract type SmoothVerticalCoordinate end

"""
    Sleve(eta_H, s)

SLEVE (Smooth LEvel VErtical) coordinate transformation for terrain-following
vertical grids in atmospheric models. Ensures that the topographic influence
decays smoothly with height, returning to flat levels well below the model top.

The transformed height is given by:
- For `η ≤ η_H`:
  `z = η·H + z_s · sinh((η_H - η) / s / η_H) / sinh(1/s)`
- For `η > η_H`:
  `z = η·H`

where `η = z_reference / H` is the normalized reference height.

# Arguments
- `eta_H`: fraction of the domain height `H` above which the
           topography influence vanishes completely.
           Must be in `(0, 1)`. Typical value: `0.7`.
- `s`: controls the rate of decay of the topographic influence with height. 
       Smaller values produce faster decay. Typical value: `0.8`.

# References
- Schär, C., Leuenberger, D., Fuhrer, O., Lüthi, D., & Frei, C. (2002)
  A new terrain-following vertical coordinate formulation for atmospheric
  prediction models
  Monthly Weather Review, 130(10), 2459–2480
  [DOI: 10.1175/1520-0493(2002)130<2459:ANTFVC>2.0.CO;2](https://doi.org/10.1175/1520-0493(2002)130<2459:ANTFVC>2.0.CO;2)
"""
struct Sleve{RealT} <: SmoothVerticalCoordinate
    eta_H::RealT  # fraction of H above which topography influence vanishes
    s::RealT           # rate of decay of the topography influence with height
end

@inline function (adapt_sleve::Sleve)(z_reference, z_topography, H)
    @unpack s, eta_H = adapt_sleve
    eta = z_reference / H
    if eta <= eta_H
        z = eta * H +
            z_topography * sinh((eta_H - eta) / s / eta_H) / sinh(one(z_reference) / s)
    else
        z = eta * H
    end
    return z
end

"""
    GalChen()

Gal-Chen & Somerville terrain-following vertical coordinate transformation.
The topographic influence decays linearly with height, reaching zero at the
top of the domain `H`.

The transformed height is given by:
  `z = z_reference + (H - z_reference) · z_s / H`

which ensures:
- At the surface (`z_reference = 0`): `z = z_s` (grid follows topography exactly)
- At the top (`z_reference = H`): `z = H` (grid is flat, independent of topography)

# References
- Gal-Chen, T. & Somerville, R. C. J. (1975)
  On the use of a coordinate transformation for the solution of the
  Navier-Stokes equations
  Journal of Computational Physics, 17(2), 209–228
  [DOI: 10.1016/0021-9991(75)90037-6](https://doi.org/10.1016/0021-9991(75)90037-6)
"""
struct GalChen <: SmoothVerticalCoordinate end

@inline function (adapt_galchen::GalChen)(z_reference, z_topography, H)
    z = z_reference + (H - z_reference) * z_topography / H
    return z
end

"""
    P4estMeshCubedSphereTopography(trees_per_face_dimension, layers, inner_radius, thickness;
                         polydeg, RealT=Float64,
                         initial_refinement_level=0, unsaved_changes=true,
                         p4est_partition_allow_for_coarsening=true, initial_topography, adapt_vertical_grid)

Build a "Cubed Sphere" mesh as `P4estMesh` with
`6 * trees_per_face_dimension^2 * layers` trees and a given topography.

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
- `initial_topography`: initial topography as a function of x, y, z coordinates, that modifies the surface spherical layer.
- `adaptive_vertical_grid`: smoothing of the vertical element size in the radial direction, to gradually restore the spherical shape. Two options are available: GalChen() or Sleve(etaH, s).        
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

# Connectivity of a cubed-sphere where we number along a column first!
function connectivity_cubed_sphere(trees_per_face_dimension, layers)
    n_cells_x = n_cells_y = trees_per_face_dimension
    n_cells_z = layers

    linear_indices = LinearIndices((layers, trees_per_face_dimension,
                                    trees_per_face_dimension, 6))

    # Vertices represent the coordinates of the forest. This is used by `p4est`
    # to write VTK files.
    # Trixi.jl doesn't use the coordinates from `p4est`, so the vertices can be empty.
    n_vertices = 0
    n_trees = 6 * n_cells_x * n_cells_y * n_cells_z
    # No edge connectivity is needed
    n_edges = 0
    # No corner connectivity is needed
    n_corners = 0
    vertices = C_NULL
    tree_to_vertex = C_NULL

    tree_to_tree = Array{p4est_topidx_t, 2}(undef, 6, n_trees)
    tree_to_face = Array{Int8, 2}(undef, 6, n_trees)

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
        for cell_z in 1:n_cells_z, cell_y in 1:n_cells_y, cell_x in 1:n_cells_x
            tree = linear_indices[cell_z, cell_x, cell_y, direction]

            # Subtract 1 because `p4est` uses zero-based indexing
            # Negative x-direction
            if cell_x > 1 # Connect to tree at the same face
                tree_to_tree[1, tree] = linear_indices[cell_z, cell_x - 1, cell_y,
                                                       direction] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 1 # This is the -x face
                target = 4
                tree_to_tree[1, tree] = linear_indices[cell_z, end, cell_y, target] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 2 # This is the +x face
                target = 3
                tree_to_tree[1, tree] = linear_indices[cell_z, end, cell_y, target] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 3 # This is the -y face
                target = 1
                tree_to_tree[1, tree] = linear_indices[cell_z, end, cell_y, target] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 4 # This is the +y face
                target = 2
                tree_to_tree[1, tree] = linear_indices[cell_z, end, cell_y, target] - 1
                tree_to_face[1, tree] = 1
            elseif direction == 5 # This is the -z face
                target = 2
                tree_to_tree[1, tree] = linear_indices[cell_z, cell_y, 1, target] - 1
                tree_to_face[1, tree] = 2
            else # direction == 6, this is the +z face
                target = 1
                tree_to_tree[1, tree] = linear_indices[cell_z, end - cell_y + 1, end,
                                                       target] - 1
                tree_to_face[1, tree] = 9 # first face dimensions are oppositely oriented, add 6
            end

            # Positive x-direction
            if cell_x < n_cells_x # Connect to tree at the same face
                tree_to_tree[2, tree] = linear_indices[cell_z, cell_x + 1, cell_y,
                                                       direction] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 1 # This is the -x face
                target = 3
                tree_to_tree[2, tree] = linear_indices[cell_z, 1, cell_y, target] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 2 # This is the +x face
                target = 4
                tree_to_tree[2, tree] = linear_indices[cell_z, 1, cell_y, target] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 3 # This is the -y face
                target = 2
                tree_to_tree[2, tree] = linear_indices[cell_z, 1, cell_y, target] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 4 # This is the +y face
                target = 1
                tree_to_tree[2, tree] = linear_indices[cell_z, 1, cell_y, target] - 1
                tree_to_face[2, tree] = 0
            elseif direction == 5 # This is the -z face
                target = 1
                tree_to_tree[2, tree] = linear_indices[cell_z, end - cell_y + 1, 1,
                                                       target] - 1
                tree_to_face[2, tree] = 8 # first face dimensions are oppositely oriented, add 6
            else # direction == 6, this is the +z face
                target = 2
                tree_to_tree[2, tree] = linear_indices[cell_z, cell_y, end, target] - 1
                tree_to_face[2, tree] = 3
            end

            # Negative y-direction
            if cell_y > 1 # Connect to tree at the same face
                tree_to_tree[3, tree] = linear_indices[cell_z, cell_x, cell_y - 1,
                                                       direction] - 1
                tree_to_face[3, tree] = 3
            elseif direction == 1
                target = 5
                tree_to_tree[3, tree] = linear_indices[cell_z, end, end - cell_x + 1,
                                                       target] - 1
                tree_to_face[3, tree] = 7 # first face dimensions are oppositely oriented, add 6
            elseif direction == 2
                target = 5
                tree_to_tree[3, tree] = linear_indices[cell_z, 1, cell_x, target] - 1
                tree_to_face[3, tree] = 0
            elseif direction == 3
                target = 5
                tree_to_tree[3, tree] = linear_indices[cell_z, end - cell_x + 1, 1,
                                                       target] - 1
                tree_to_face[3, tree] = 8 # first face dimensions are oppositely oriented, add 6
            elseif direction == 4
                target = 5
                tree_to_tree[3, tree] = linear_indices[cell_z, cell_x, end, target] - 1
                tree_to_face[3, tree] = 3
            elseif direction == 5
                target = 3
                tree_to_tree[3, tree] = linear_indices[cell_z, end - cell_x + 1, 1,
                                                       target] - 1
                tree_to_face[3, tree] = 8 # first face dimensions are oppositely oriented, add 6
            else # direction == 6
                target = 3
                tree_to_tree[3, tree] = linear_indices[cell_z, cell_x, end, target] - 1
                tree_to_face[3, tree] = 3
            end

            # Positive y-direction
            if cell_y < n_cells_y # Connect to tree at the same face
                tree_to_tree[4, tree] = linear_indices[cell_z, cell_x, cell_y + 1,
                                                       direction] - 1
                tree_to_face[4, tree] = 2
            elseif direction == 1
                target = 6
                tree_to_tree[4, tree] = linear_indices[cell_z, 1, end - cell_x + 1,
                                                       target] - 1
                tree_to_face[4, tree] = 6 # first face dimensions are oppositely oriented, add 6
            elseif direction == 2
                target = 6
                tree_to_tree[4, tree] = linear_indices[cell_z, end, cell_x, target] - 1
                tree_to_face[4, tree] = 1
            elseif direction == 3
                target = 6
                tree_to_tree[4, tree] = linear_indices[cell_z, cell_x, 1, target] - 1
                tree_to_face[4, tree] = 2
            elseif direction == 4
                target = 6
                tree_to_tree[4, tree] = linear_indices[cell_z, end - cell_x + 1, end,
                                                       target] - 1
                tree_to_face[4, tree] = 9 # first face dimensions are oppositely oriented, add 6
            elseif direction == 5
                target = 4
                tree_to_tree[4, tree] = linear_indices[cell_z, cell_x, 1, target] - 1
                tree_to_face[4, tree] = 2
            else # direction == 6
                target = 4
                tree_to_tree[4, tree] = linear_indices[cell_z, end - cell_x + 1, end,
                                                       target] - 1
                tree_to_face[4, tree] = 9 # first face dimensions are oppositely oriented, add 6
            end

            # Negative z-direction
            if cell_z > 1
                tree_to_tree[5, tree] = linear_indices[cell_z - 1, cell_x, cell_y,
                                                       direction] - 1
                tree_to_face[5, tree] = 5
            else # Non-periodic boundary, tree and face point to themselves (zero-based indexing)
                tree_to_tree[5, tree] = tree - 1
                tree_to_face[5, tree] = 4
            end

            # Positive z-direction
            if cell_z < n_cells_z
                tree_to_tree[6, tree] = linear_indices[cell_z + 1, cell_x, cell_y,
                                                       direction] - 1
                tree_to_face[6, tree] = 4
            else # Non-periodic boundary, tree and face point to themselves (zero-based indexing)
                tree_to_tree[6, tree] = tree - 1
                tree_to_face[6, tree] = 5
            end
        end
    end

    tree_to_edge = C_NULL
    # `p4est` docs: "in trivial cases it is just a pointer to a p4est_topix value of 0."
    # We don't need edge connectivity, so this is a trivial case.
    ett_offset = zeros(p4est_topidx_t, 1)
    edge_to_tree = C_NULL
    edge_to_edge = C_NULL

    tree_to_corner = C_NULL
    # `p4est` docs: "in trivial cases it is just a pointer to a p4est_topix value of 0."
    # We don't need corner connectivity, so this is a trivial case.
    ctt_offset = zeros(p4est_topidx_t, 1)

    corner_to_tree = C_NULL
    corner_to_corner = C_NULL

    connectivity = p8est_connectivity_new_copy(n_vertices, n_trees, n_corners, n_edges,
                                               vertices, tree_to_vertex,
                                               tree_to_tree, tree_to_face,
                                               tree_to_edge, ett_offset,
                                               edge_to_tree, edge_to_edge,
                                               tree_to_corner, ctt_offset,
                                               corner_to_tree, corner_to_corner)

    @assert p8est_connectivity_is_valid(connectivity) == 1

    return connectivity
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
    lat = asin(unit_pt[3] / r)
    lon = atan(unit_pt[2], unit_pt[1])
    z_topography = initial_topography(lat, lon)

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

    linear_indices = LinearIndices((n_cells_z, n_cells_x, n_cells_y, 6))

    dx = 2 / n_cells_x
    dy = 2 / n_cells_y
    dz = 2 / n_cells_z

    for direction in 1:6
        for cell_z in 1:n_cells_z, cell_y in 1:n_cells_y, cell_x in 1:n_cells_x
            tree = linear_indices[cell_z, cell_x, cell_y, direction]

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
