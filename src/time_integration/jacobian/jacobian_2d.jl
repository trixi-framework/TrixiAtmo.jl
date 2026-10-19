# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Vertical layout of the degrees of freedom and metric terms of the vertical Jacobians of
# the `CompressibleEulerPotentialTemperatureEquationsWithGravity2D` on a `P4estMesh{2}`,
# computed once on the CPU

# The trees of a `P4estMesh{2}` created with `trees_per_dimension` are ordered with the
# first (horizontal) coordinate fastest, such that the layers are the rows of trees in the
# direction of the second (vertical) coordinate. The number of layers follows from the
# connectivity of the trees.
function number_of_layers(mesh::Trixi.P4estMesh{2})
    # The elements of a column are not consecutive in the order of the trees, so they would
    # be distributed over several ranks
    if mesh isa Trixi.P4estMeshParallel{2}
        throw(ArgumentError("the vertical Jacobians do not support parallel " *
                            "`P4estMesh{2}`"))
    end
    connectivity = mesh.p4est.connectivity
    ntrees_ = Int(connectivity.num_trees[])
    # Neighbor of the zero-based `tree` at its zero-based `face`, where `3 => +y`
    neighbor(tree, face) = Int(connectivity.tree_to_tree[4 * tree + face + 1])

    # Walk upwards from the first tree to the upper boundary, where the neighbor of a
    # tree is the tree itself
    nlayers = 1
    tree = 0
    while neighbor(tree, 3) != tree
        tree = neighbor(tree, 3)
        if tree == 0
            throw(ArgumentError("the vertical Jacobians require a mesh that is not " *
                                "periodic in the vertical direction"))
        end
        nlayers += 1
    end

    nhorizontal = ntrees_ ÷ nlayers
    for tree in 0:(ntrees_ - nhorizontal - 1)
        if neighbor(tree, 3) != tree + nhorizontal
            throw(ArgumentError("the vertical Jacobians require a `P4estMesh{2}` " *
                                "created with `trees_per_dimension`"))
        end
    end
    return nlayers
end

# In 2D, the layers are the rows of trees (see `number_of_layers`), with the lowest layer
# first
function compute_vertical_mapping(::Trixi.P4estMesh{2}, nnodes_, nelements_, nlayers)
    vertical_line_indices = zeros(Int, nnodes_, nnodes_, nelements_)
    vertical_element_indices = zeros(Int, nnodes_, nnodes_, nelements_)

    nhorizontal = nelements_ ÷ nlayers
    for element in 1:nelements_
        layer = (element - 1) ÷ nhorizontal + 1
        column_offset = ((element - 1) % nhorizontal) * nnodes_
        for j in 1:nnodes_, i in 1:nnodes_
            vertical_line_indices[i, j, element] = column_offset + i
            vertical_element_indices[i, j, element] = layer
        end
    end

    return vertical_line_indices, vertical_element_indices
end

# In 2D, a column collects `nnodes` nodes of each element instead of `nnodes^2`
function compute_element_index_maps(vertical_line_indices::AbstractArray{<:Any, 3},
                                    vertical_element_indices, nnodes_)
    nelements_ = size(vertical_line_indices, 3)
    horizontal_index_of_element = [(vertical_line_indices[1, 1, element] - 1) ÷
                                   nnodes_ + 1 for element in 1:nelements_]
    vertical_index_of_element = [vertical_element_indices[1, 1, element]
                                 for element in 1:nelements_]
    for element in 1:nelements_, j in 1:nnodes_, i in 1:nnodes_
        @assert vertical_line_indices[i, j, element] ==
                (horizontal_index_of_element[element] - 1) * nnodes_ + i
        @assert vertical_element_indices[i, j, element] ==
                vertical_index_of_element[element]
    end
    return horizontal_index_of_element, vertical_index_of_element
end

# In 2D, `vertical_frame[1:2, ...]` holds the unit contravariant direction `nu = Ja^2 / n`
# of the second reference coordinate, oriented upwards (against gravity in the negative
# y-direction), and `vertical_frame[3, ...]` the normalized norm `n = |Ja^2|`
function compute_vertical_frame(::Trixi.P4estMesh{2}, semi)
    (; solver) = semi
    contravariant_vectors = Adapt.adapt(Array,
                                        semi.cache.elements.contravariant_vectors)
    vertical_frame = zeros(eltype(contravariant_vectors), 3,
                           size(contravariant_vectors)[3:end]...)
    area_max = zero(eltype(contravariant_vectors))
    for element in axes(contravariant_vectors, 5)
        for j in eachnode(solver), i in eachnode(solver)
            Ja2 = get_contravariant_vector(2, contravariant_vectors, i, j, element)
            area = norm(Ja2)
            direction = sign(Ja2[2]) * Ja2 / area
            vertical_frame[1:2, i, j, element] .= direction
            vertical_frame[3, i, j, element] = area
            area_max = max(area_max, area)
        end
    end
    vertical_frame[3, :, :, :] ./= area_max
    return vertical_frame
end

# In 2D, the metric factor is `n / J` with the norm `n = |Ja^2|`
function compute_metric_factor(::Trixi.P4estMesh{2}, ncolumns, nnodes_, nlayers, semi,
                               vertical_line_indices,
                               vertical_element_indices)
    contravariant_vectors = Adapt.adapt(Array,
                                        semi.cache.elements.contravariant_vectors)
    inverse_jacobian = Adapt.adapt(Array, semi.cache.elements.inverse_jacobian)
    metric_factor = zeros(eltype(inverse_jacobian), ncolumns, nnodes_, nlayers)
    for element in axes(inverse_jacobian, 3), j in 1:nnodes_, i in 1:nnodes_
        Ja2 = get_contravariant_vector(2, contravariant_vectors, i, j, element)
        column = vertical_line_indices[i, j, element]
        layer = vertical_element_indices[i, j, element]
        metric_factor[column, j, layer] = norm(Ja2) * inverse_jacobian[i, j, element]
    end
    return metric_factor
end

# In 2D, gravity acts in the negative y-direction, i.e., the pointwise gravity is
# `g * (nu . y_hat)`. It reduces to `g` if the second reference coordinate is aligned with
# `y`.
function compute_gravity(::Trixi.P4estMesh{2}, ncolumns, nnodes_, nlayers, semi,
                         vertical_line_indices, vertical_element_indices)
    (; equations) = semi
    contravariant_vectors = Adapt.adapt(Array,
                                        semi.cache.elements.contravariant_vectors)
    gravity = zeros(eltype(contravariant_vectors), ncolumns, nnodes_, nlayers)
    for element in axes(contravariant_vectors, 5), j in 1:nnodes_, i in 1:nnodes_
        Ja2 = get_contravariant_vector(2, contravariant_vectors, i, j, element)
        # |nu . y_hat| with the upward orientation of `compute_vertical_frame`
        projection = abs(Ja2[2]) / norm(Ja2)
        column = vertical_line_indices[i, j, element]
        layer = vertical_element_indices[i, j, element]
        gravity[column, j, layer] = equations.g * projection
    end
    return gravity
end
end # @muladd
