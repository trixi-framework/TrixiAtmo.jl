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

# Column and layer of every element, used to copy data between Trixi.jl's layout and the
# columns. The formula is specific to the ordering of the trees of the
# `P4estMeshCubedSphereTopography` (see `calc_tree_node_coordinates!` in
# `src/meshes/p4est_cubed_sphere_topography.jl`), where the layers of a column are
# consecutive, with the lowest layer first; other meshes are rejected by `number_of_layers`,
# and the ordering is checked by `check_vertical_ordering`.
function compute_element_index_maps(nlayers, solver, cache)
    horizontal_index_of_element = [(element - 1) ÷ nlayers + 1
                                   for element in eachelement(solver, cache)]
    vertical_index_of_element = [(element - 1) % nlayers + 1
                                 for element in eachelement(solver, cache)]
    return horizontal_index_of_element, vertical_index_of_element
end

# Index of the column of the nodes `(i, j, k)` of `element` (for all `k`), i.e., the first
# index of the arrays laid out `[column, ...]`
Base.@propagate_inbounds function column_index(horizontal_index_of_element, nnodes_,
                                               i, j, element)
    return (horizontal_index_of_element[element] - 1) * nnodes_ * nnodes_ +
           (j - 1) * nnodes_ + i
end

# Check the assumptions of the layout of the columns on the geometry of the mesh: the third
# reference coordinate points outwards (`Ja^3 . x > 0`), and the last nodes (`k = nnodes`)
# of every element except the top ones coincide with the first nodes (`k = 1`) of the next
# element, i.e., the elements of a column are consecutive with the lowest layer first
function check_vertical_ordering(nlayers, semi)
    (; solver, equations, cache) = semi
    nnodes_ = nnodes(solver)
    nelements_ = nelements(solver, cache)
    if nelements_ % nlayers != 0
        throw(ArgumentError("the number of elements $nelements_ is not a multiple of " *
                            "the number of layers $nlayers"))
    end
    node_coordinates = Adapt.adapt(Array, cache.elements.node_coordinates)
    contravariant_vectors = Adapt.adapt(Array, cache.elements.contravariant_vectors)
    RealT = eltype(node_coordinates)
    tolerance = sqrt(eps(RealT))
    for element in eachelement(solver, cache)
        for k in eachnode(solver), j in eachnode(solver), i in eachnode(solver)
            Ja3 = get_contravariant_vector(3, contravariant_vectors, i, j, k, element)
            x = get_node_coords(node_coordinates, equations, solver, i, j, k, element)
            if !(dot(Ja3, x) > 0)
                throw(ArgumentError("the third reference coordinate of element " *
                                    "$element must point outwards"))
            end
        end
        # The top element of a column has no element above it in the same column
        (element - 1) % nlayers + 1 == nlayers && continue
        for j in eachnode(solver), i in eachnode(solver)
            x_top = get_node_coords(node_coordinates, equations, solver, i, j, nnodes_,
                                    element)
            x_bottom = get_node_coords(node_coordinates, equations, solver, i, j, 1,
                                       element + 1)
            if !isapprox(x_top, x_bottom; rtol = tolerance)
                throw(ArgumentError("the elements $element and $(element + 1) are not " *
                                    "consecutive layers of a column"))
            end
        end
    end
    return nothing
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
    norm_Ja3_max = zero(eltype(node_coordinates))
    for element in eachelement(solver, cache)
        for k in eachnode(solver), j in eachnode(solver), i in eachnode(solver)
            Ja3 = get_contravariant_vector(3, contravariant_vectors, i, j, k, element)
            norm_Ja3 = norm(Ja3)
            x = get_node_coords(node_coordinates, equations, solver, i, j, k, element)
            # `Ja3` points outwards, see `check_vertical_ordering`
            vertical_frame[1:3, i, j, k, element] .= Ja3 / norm_Ja3
            vertical_frame[4, i, j, k, element] = norm_Ja3
            norm_Ja3_max = max(norm_Ja3_max, norm_Ja3)
        end
    end
    vertical_frame[4, :, :, :, :] ./= norm_Ja3_max
    return vertical_frame
end

# Metric factor `n / J` of the vertical operator of every node, with the norm
# `n = |Ja^3|` of the contravariant vector of the third reference coordinate and the
# Jacobian `J`, laid out `[column, node, layer]`
function compute_metric_factor(ncolumns, nlayers, semi, horizontal_index_of_element,
                               vertical_index_of_element)
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
        column = column_index(horizontal_index_of_element, nnodes_, i, j, element)
        layer = vertical_index_of_element[element]
        metric_factor[column, k, layer] = norm(Ja3) *
                                          inverse_jacobian[i, j, k, element]
    end
    return metric_factor
end

# Consistent pointwise gravity `g * (nu . r_hat)`, the projection of the gravitational
# acceleration on the unit contravariant direction `nu` of the third reference coordinate
# (see `compute_vertical_frame`).
function compute_gravity(ncolumns, nlayers, semi, horizontal_index_of_element,
                         vertical_index_of_element)
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
            # `nu . r_hat` is positive since `Ja3` points outwards, see
            # `check_vertical_ordering`
            projection = dot(Ja3, x) / (norm(Ja3) * norm(x))
            column = column_index(horizontal_index_of_element, nnodes_, i, j, element)
            layer = vertical_index_of_element[element]
            gravity[column, k, layer] = equations.g * projection
        end
    end
    return gravity
end
end # @muladd
