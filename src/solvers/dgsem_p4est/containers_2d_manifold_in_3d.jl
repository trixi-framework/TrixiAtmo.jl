@muladd begin
#! format: noindent

# New p4est element container that allows the use of a PtrArray for the contravariant_vectors
mutable struct P4estElementContainerPtrArray{NDIMS, RealT <: Real, uEltype <: Real,
                                             NDIMSP1,
                                             NDIMSP2, NDIMSP3,
                                             ContravariantVectors <:
                                             AbstractArray{RealT, NDIMSP3}} <:
               Trixi.AbstractContainer
    # Physical coordinates at each node
    node_coordinates::Array{RealT, NDIMSP2}   # [orientation, node_i, node_j, node_k, element]
    # Jacobian matrix of the transformation
    # [jacobian_i, jacobian_j, node_i, node_j, node_k, element] where jacobian_i is the first index of the Jacobian matrix,...
    jacobian_matrix::Array{RealT, NDIMSP3}
    # Contravariant vectors, scaled by J, in Kopriva's blue book called Ja^i_n (i index, n dimension)
    contravariant_vectors::ContravariantVectors  # [dimension, index, node_i, node_j, node_k, element]
    # 1/J where J is the Jacobian determinant (determinant of Jacobian matrix)
    inverse_jacobian::Array{RealT, NDIMSP1}   # [node_i, node_j, node_k, element]
    # Buffer for calculated surface flux
    surface_flux_values::Array{uEltype, NDIMSP2} # [variable, i, j, direction, element]

    # internal `resize!`able storage
    _node_coordinates::Vector{RealT}
    _jacobian_matrix::Vector{RealT}
    _contravariant_vectors::Vector{RealT}
    _inverse_jacobian::Vector{RealT}
    _surface_flux_values::Vector{uEltype}
end

@inline function Trixi.nelements(elements::P4estElementContainerPtrArray)
    size(elements.node_coordinates, ndims(elements) + 2)
end
@inline Base.ndims(::P4estElementContainerPtrArray{NDIMS}) where {NDIMS} = NDIMS
@inline function Base.eltype(::P4estElementContainerPtrArray{NDIMS, RealT, uEltype}) where {
                                                                                            NDIMS,
                                                                                            RealT,
                                                                                            uEltype
                                                                                            }
    uEltype
end

# Create element container and initialize element data.
# This function dispatches on the dimensions of the mesh and the equation (AbstractEquations{3})
function Trixi.init_elements(mesh::Union{P4estMesh{2, 3, RealT},
                                         T8codeMesh{2}},
                             equations::AbstractEquations{3},
                             basis,
                             ::Type{uEltype}) where {RealT <: Real, uEltype <: Real}
    nelements = Trixi.ncells(mesh)

    NDIMS = 2 # dimension of the manifold
    NDIMS_AMBIENT = 3 # dimension of the ambient space

    _node_coordinates = Vector{RealT}(undef,
                                      NDIMS_AMBIENT * nnodes(basis)^NDIMS * nelements)
    node_coordinates = Trixi.unsafe_wrap(Array, pointer(_node_coordinates),
                                         (NDIMS_AMBIENT,
                                          ntuple(_ -> nnodes(basis), NDIMS)...,
                                          nelements))

    _jacobian_matrix = Vector{RealT}(undef,
                                     NDIMS_AMBIENT * NDIMS * nnodes(basis)^NDIMS *
                                     nelements)
    jacobian_matrix = Trixi.unsafe_wrap(Array, pointer(_jacobian_matrix),
                                        (NDIMS_AMBIENT, NDIMS,
                                         ntuple(_ -> nnodes(basis), NDIMS)...,
                                         nelements))

    _contravariant_vectors = Vector{RealT}(undef,
                                           NDIMS_AMBIENT^2 * nnodes(basis)^NDIMS *
                                           nelements)
    contravariant_vectors = PtrArray(pointer(_contravariant_vectors),
                                     (Trixi.StaticInt(NDIMS_AMBIENT), NDIMS_AMBIENT,
                                      ntuple(_ -> nnodes(basis), NDIMS)...,
                                      nelements))

    _inverse_jacobian = Vector{RealT}(undef, nnodes(basis)^NDIMS * nelements)
    inverse_jacobian = Trixi.unsafe_wrap(Array, pointer(_inverse_jacobian),
                                         (ntuple(_ -> nnodes(basis), NDIMS)...,
                                          nelements))

    _surface_flux_values = Vector{uEltype}(undef,
                                           nvariables(equations) *
                                           nnodes(basis)^(NDIMS - 1) * (NDIMS * 2) *
                                           nelements)
    surface_flux_values = Trixi.unsafe_wrap(Array, pointer(_surface_flux_values),
                                            (nvariables(equations),
                                             ntuple(_ -> nnodes(basis), NDIMS - 1)...,
                                             NDIMS * 2, nelements))

    elements = P4estElementContainerPtrArray{NDIMS, RealT, uEltype, NDIMS + 1,
                                             NDIMS + 2,
                                             NDIMS + 3, typeof(contravariant_vectors)}(node_coordinates,
                                                                                       jacobian_matrix,
                                                                                       contravariant_vectors,
                                                                                       inverse_jacobian,
                                                                                       surface_flux_values,
                                                                                       _node_coordinates,
                                                                                       _jacobian_matrix,
                                                                                       _contravariant_vectors,
                                                                                       _inverse_jacobian,
                                                                                       _surface_flux_values)

    init_elements_2d_manifold_in_3d!(elements, mesh, basis)

    return elements
end

function init_elements_2d_manifold_in_3d!(elements,
                                          mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                          basis::LobattoLegendreBasis)
    (; node_coordinates, jacobian_matrix,
    contravariant_vectors, inverse_jacobian) = elements

    calc_node_coordinates_2d_shell!(node_coordinates, mesh, basis)

    for element in 1:Trixi.ncells(mesh)

        # Compute Jacobian matrix as usual
        Trixi.calc_jacobian_matrix!(jacobian_matrix, element, node_coordinates, basis)

        # Compute contravariant vectors
        calc_contravariant_vectors_2d_shell!(contravariant_vectors,
                                             element,
                                             jacobian_matrix,
                                             node_coordinates,
                                             basis)

        # Compute the inverse Jacobian as the norm of the cross product of the covariant vectors
        for j in eachnode(basis), i in eachnode(basis)
            inverse_jacobian[i, j, element] = 1 /
                                              norm(Trixi.cross(jacobian_matrix[:, 1, i,
                                                                               j,
                                                                               element],
                                                               jacobian_matrix[:, 2, i,
                                                                               j,
                                                                               element]))
        end
    end

    return nothing
end

# Interpolate tree_node_coordinates to each quadrant at the nodes of the specified basis
function calc_node_coordinates_2d_shell!(node_coordinates,
                                         mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                         basis::LobattoLegendreBasis)
    # Hanging nodes will cause holes in the mesh if its polydeg is higher
    # than the polydeg of the solver.
    @assert length(basis.nodes)>=length(mesh.nodes) "The solver can't have a lower polydeg than the mesh"
    @assert size(mesh.tree_node_coordinates, 1)==3 "Shell must use 3D node coordinates"
    calc_node_coordinates_2d_shell!(node_coordinates, mesh, basis.nodes)
end

# Interpolate tree_node_coordinates to each quadrant at the specified nodes
function calc_node_coordinates_2d_shell!(node_coordinates,
                                         mesh::P4estMesh{2},
                                         nodes::AbstractVector)
    # We use `StrideArray`s here since these buffers are used in performance-critical
    # places and the additional information passed to the compiler makes them faster
    # than native `Array`s.
    tmp1 = Trixi.StrideArray(undef, real(mesh),
                             Trixi.StaticInt(3),
                             Trixi.static_length(nodes),
                             Trixi.static_length(mesh.nodes))
    matrix1 = Trixi.StrideArray(undef, real(mesh),
                                Trixi.static_length(nodes),
                                Trixi.static_length(mesh.nodes))
    matrix2 = similar(matrix1)
    baryweights_in = Trixi.barycentric_weights(mesh.nodes)

    # Macros from `p4est`
    p4est_root_len = 1 << Trixi.P4EST_MAXLEVEL
    p4est_quadrant_len(l) = 1 << (Trixi.P4EST_MAXLEVEL - l)

    trees = Trixi.unsafe_wrap_sc(Trixi.p4est_tree_t, mesh.p4est.trees)

    for tree in eachindex(trees)
        offset = trees[tree].quadrants_offset
        quadrants = Trixi.unsafe_wrap_sc(Trixi.p4est_quadrant_t, trees[tree].quadrants)

        for i in eachindex(quadrants)
            element = offset + i
            quad = quadrants[i]

            quad_length = p4est_quadrant_len(quad.level) / p4est_root_len

            nodes_out_x = 2 * (quad_length * 1 / 2 * (nodes .+ 1) .+
                           quad.x / p4est_root_len) .- 1
            nodes_out_y = 2 * (quad_length * 1 / 2 * (nodes .+ 1) .+
                           quad.y / p4est_root_len) .- 1
            Trixi.polynomial_interpolation_matrix!(matrix1, mesh.nodes, nodes_out_x,
                                                   baryweights_in)
            Trixi.polynomial_interpolation_matrix!(matrix2, mesh.nodes, nodes_out_y,
                                                   baryweights_in)

            Trixi.multiply_dimensionwise!(view(node_coordinates, :, :, :, element),
                                          matrix1, matrix2,
                                          view(mesh.tree_node_coordinates, :, :, :,
                                               tree),
                                          tmp1)
        end
    end

    return node_coordinates
end

# This only works for a sphere
function calc_contravariant_vectors_2d_shell!(contravariant_vectors::AbstractArray{<:Any,
                                                                                   5},
                                              element,
                                              jacobian_matrix, node_coordinates,
                                              basis::LobattoLegendreBasis)
    @unpack derivative_matrix = basis

    for j in eachnode(basis), i in eachnode(basis)
        radius = sqrt(node_coordinates[1, i, j, element]^2 +
                      node_coordinates[2, i, j, element]^2 +
                      node_coordinates[3, i, j, element]^2)

        for n in 1:3
            # (n, m, l) cyclic
            m = (n % 3) + 1
            l = ((n + 1) % 3) + 1

            contravariant_vectors[n, 1, i, j, element] = (jacobian_matrix[m, 2, i, j,
                                                                          element] *
                                                          node_coordinates[l, i, j,
                                                                           element]
                                                          -
                                                          jacobian_matrix[l, 2, i, j,
                                                                          element] *
                                                          node_coordinates[m, i, j,
                                                                           element]) /
                                                         radius

            contravariant_vectors[n, 2, i, j, element] = (jacobian_matrix[l, 1, i, j,
                                                                          element] *
                                                          node_coordinates[m, i, j,
                                                                           element]
                                                          -
                                                          jacobian_matrix[m, 1, i, j,
                                                                          element] *
                                                          node_coordinates[l, i, j,
                                                                           element]) /
                                                         radius

            contravariant_vectors[n, 3, i, j, element] = (jacobian_matrix[m, 1, i, j,
                                                                          element] *
                                                          jacobian_matrix[l, 2, i, j,
                                                                          element]
                                                          -
                                                          jacobian_matrix[m, 2, i, j,
                                                                          element] *
                                                          jacobian_matrix[l, 1, i, j,
                                                                          element]) /
                                                         radius
        end
    end

    return contravariant_vectors
end
end # @muladd
