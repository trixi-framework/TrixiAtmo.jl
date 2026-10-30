# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Vertical layout of the degrees of freedom and metric terms of the vertical Jacobians of
# the `CompressibleEulerPotentialTemperatureEquationsWithGravity3D` on a
# `P4estMeshCubedSphereTopography`, computed once on the CPU

# The number of layers of elements of a `P4estMeshCubedSphereTopography` is stored in the
# registry of the mesh (see `src/meshes/p4est_cubed_sphere_topography.jl`). For parallel
# meshes, the columns must be kept together on one rank.
function number_of_layers(mesh::Trixi.P4estMesh{3})
    key = UInt64(Trixi.P4est.PointerWrappers.pointer(mesh.p4est))
    if !haskey(MeshContextRegistry.REGISTRY, key)
        throw(ArgumentError("the vertical Jacobians require a " *
                            "`P4estMeshCubedSphereTopography`"))
    end
    ctx = MeshContextRegistry.REGISTRY[key]::ColumnPartitionContext
    if mesh isa Trixi.P4estMeshParallel{3} && !ctx.keep_columns_together
        throw(ArgumentError("the vertical Jacobians require " *
                            "`keep_columns_together = true` for parallel meshes"))
    end
    return Int(ctx.num_layers_per_column)
end

# Column and layer of every node. The trees of the cubed sphere are ordered such that the
# layers of a column are consecutive, with the lowest layer first.
function compute_vertical_mapping(nlayers, solver, cache)
    nnodes_ = nnodes(solver)
    nelements_ = nelements(solver, cache)
    vertical_line_indices = zeros(Int, nnodes_, nnodes_, nnodes_, nelements_)
    vertical_element_indices = zeros(Int, nnodes_, nnodes_, nnodes_, nelements_)

    for element in eachelement(solver, cache)
        layer = (element - 1) % nlayers + 1
        column_offset = ((element - 1) ÷ nlayers) * nnodes_^2
        for k in eachnode(solver), j in eachnode(solver), i in eachnode(solver)
            vertical_line_indices[i, j, k, element] = column_offset +
                                                      (j - 1) * nnodes_ +
                                                      i
            vertical_element_indices[i, j, k, element] = layer
        end
    end

    return vertical_line_indices, vertical_element_indices
end

# Per-element maps used to copy data between Trixi.jl's layout and the columns: the
# horizontal column index and the layer, one per element instead of one per node.
function compute_element_index_maps(vertical_line_indices, vertical_element_indices,
                                    solver, cache)
    nnodes_ = nnodes(solver)
    horizontal_index_of_element = [(vertical_line_indices[1, 1, 1, element] - 1) ÷
                                   nnodes_^2 + 1
                                   for element in eachelement(solver, cache)]
    vertical_index_of_element = [vertical_element_indices[1, 1, 1, element]
                                 for element in eachelement(solver, cache)]
    for element in eachelement(solver, cache), k in eachnode(solver),
        j in eachnode(solver), i in eachnode(solver)

        @assert vertical_line_indices[i, j, k, element] ==
                (horizontal_index_of_element[element] - 1) * nnodes_^2 +
                (j - 1) * nnodes_ + i
        @assert vertical_element_indices[i, j, k, element] ==
                vertical_index_of_element[element]
    end
    return horizontal_index_of_element, vertical_index_of_element
end

# Vertical frame of every node: `vertical_frame[1:3, ...]` holds the unit contravariant
# direction `nu = Ja^3 / n` of the third reference coordinate, oriented outwards, and
# `vertical_frame[4, ...]` the norm `n = |Ja^3|` normalized by its global maximum (a
# global constant cancels in the similarity transformation `N` of the copies in
# `jacobian_3d_gpu.jl`)
function compute_vertical_frame(semi)
    (; solver, equations, cache) = semi
    node_coordinates = Adapt.adapt(Array, semi.cache.elements.node_coordinates)
    contravariant_vectors = Adapt.adapt(Array,
                                        semi.cache.elements.contravariant_vectors)
    vertical_frame = zeros(eltype(node_coordinates), 4,
                           size(node_coordinates)[2:end]...)
    area_max = zero(eltype(node_coordinates))
    for element in eachelement(solver, cache)
        for k in eachnode(solver), j in eachnode(solver), i in eachnode(solver)
            Ja3 = get_contravariant_vector(3, contravariant_vectors, i, j, k, element)
            area = norm(Ja3)
            x = get_node_coords(node_coordinates, equations, solver, i, j, k, element)
            direction = sign(dot(Ja3, x)) * Ja3 / area
            vertical_frame[1:3, i, j, k, element] .= direction
            vertical_frame[4, i, j, k, element] = area
            area_max = max(area_max, area)
        end
    end
    vertical_frame[4, :, :, :, :] ./= area_max
    return vertical_frame
end

# Metric factor `n / J` of the vertical operator of every node, with the norm
# `n = |Ja^3|` of the contravariant vector of the third reference coordinate and the
# Jacobian `J`, laid out `[column, node, layer]`
function compute_metric_factor(ncolumns, nlayers, semi, vertical_line_indices,
                               vertical_element_indices)
    (; solver, cache) = semi
    nnodes_ = nnodes(solver)
    contravariant_vectors = Adapt.adapt(Array,
                                        semi.cache.elements.contravariant_vectors)
    inverse_jacobian = Adapt.adapt(Array, semi.cache.elements.inverse_jacobian)
    metric_factor = zeros(eltype(inverse_jacobian), ncolumns, nnodes_,
                          nlayers)
    for element in eachelement(solver, cache), k in eachnode(solver),
        j in eachnode(solver), i in eachnode(solver)

        Ja3 = get_contravariant_vector(3, contravariant_vectors, i, j, k, element)
        column = vertical_line_indices[i, j, k, element]
        layer = vertical_element_indices[i, j, k, element]
        metric_factor[column, k, layer] = norm(Ja3) *
                                          inverse_jacobian[i, j, k, element]
    end
    return metric_factor
end

# Consistent pointwise gravity `g * (nu . r_hat)`, the projection of the gravitational
# acceleration on the unit contravariant direction `nu` of the third reference coordinate
# (see `compute_vertical_frame`).
function compute_gravity(ncolumns, nlayers, semi, vertical_line_indices,
                         vertical_element_indices)
    (; solver, equations, cache) = semi
    nnodes_ = nnodes(solver)
    node_coordinates = Adapt.adapt(Array, semi.cache.elements.node_coordinates)
    contravariant_vectors = Adapt.adapt(Array,
                                        semi.cache.elements.contravariant_vectors)
    gravity = zeros(eltype(node_coordinates), ncolumns, nnodes_, nlayers)
    for element in eachelement(solver, cache)
        for k in eachnode(solver), j in eachnode(solver), i in eachnode(solver)
            Ja3 = get_contravariant_vector(3, contravariant_vectors, i, j, k, element)
            x = get_node_coords(node_coordinates, equations, solver, i, j, k, element)
            # |nu . r_hat| with the outward orientation of `compute_vertical_frame`
            projection = abs(dot(Ja3, x)) / (norm(Ja3) * norm(x))
            column = vertical_line_indices[i, j, k, element]
            layer = vertical_element_indices[i, j, k, element]
            gravity[column, k, layer] = equations.g * projection
        end
    end
    return gravity
end
end # @muladd
