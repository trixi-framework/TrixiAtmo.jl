using Trixi: connectivity_cubed_sphere, new_p4est

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
