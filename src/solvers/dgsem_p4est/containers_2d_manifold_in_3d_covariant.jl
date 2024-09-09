@muladd begin
#! format: noindent

# Specialized element container for partial differential equations in covariant form
mutable struct P4estElementContainerCovariant{NDIMS, RealT <: Real, uEltype <: Real,
                                              NDIMSP1, NDIMSP2,
                                              NDIMSP3} <: Trixi.AbstractContainer
    # Physical Cartesian coordinates at each node
    node_coordinates::Array{RealT, NDIMSP2}  # [orientation, node_i, node_j, node_k, element]
    # Jacobian matrix of the transformation
    # [jacobian_i, jacobian_j, node_i, node_j, node_k, element] where jacobian_i is the first index of the Jacobian matrix,...
    jacobian_matrix::Array{RealT, NDIMSP3}
    # Transformation matrix from contravariant to zonal/meridional velocity or momentum
    transform_matrix::Array{RealT, NDIMSP3} # indexed same way as Jacobian matrix
    # Inverse of the transformation matrix
    inverse_transform_matrix::Array{RealT, NDIMSP3} # indexed same way as Jacobian matrix
    # Determinant of transformation matrix
    jacobian::Array{RealT, NDIMSP1} # [node_i, node_j, node_k, element]
    # Determinant of inverse transformation matrix
    inverse_jacobian::Array{RealT, NDIMSP1}  # [node_i, node_j, node_k, element]
    # Buffer for calculated surface flux
    surface_flux_values::Array{uEltype, NDIMSP2} # [variable, i, j, direction, element]

    # internal `resize!`able storage
    _node_coordinates::Vector{RealT}
    _jacobian_matrix::Vector{RealT}
    _transform_matrix::Vector{RealT}
    _inverse_transform_matrix::Vector{RealT}
    _jacobian::Vector{RealT}
    _inverse_jacobian::Vector{RealT}
    _surface_flux_values::Vector{uEltype}
end

@inline function Trixi.nelements(elements::P4estElementContainerCovariant)
    size(elements.node_coordinates, ndims(elements) + 2)
end
@inline Base.ndims(::P4estElementContainerCovariant{NDIMS}) where {NDIMS} = NDIMS
@inline function Base.eltype(::P4estElementContainerCovariant{NDIMS, RealT,
                                                              uEltype}) where {NDIMS,
                                                                               RealT,
                                                                               uEltype}
    return uEltype
end

# Create element container and initialize element data.
# This function dispatches on the dimensions of the mesh and the equation type 
function Trixi.init_elements(mesh::P4estMesh{NDIMS},
                             equations::AbstractCovariantEquations{NDIMS},
                             basis, ::Type{uEltype}) where {NDIMS, uEltype <: Real}
    RealT = real(mesh)
    NDIMS_AMBIENT = size(mesh.tree_node_coordinates, 1)

    nelements = Trixi.ncells(mesh)

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

    _transform_matrix = Vector{RealT}(undef,
                                      NDIMS * NDIMS * nnodes(basis)^NDIMS *
                                      nelements)
    transform_matrix = Trixi.unsafe_wrap(Array, pointer(_transform_matrix),
                                         (NDIMS, NDIMS,
                                          ntuple(_ -> nnodes(basis), NDIMS)...,
                                          nelements))

    _inverse_transform_matrix = Vector{RealT}(undef,
                                              NDIMS * NDIMS * nnodes(basis)^NDIMS *
                                              nelements)
    inverse_transform_matrix = Trixi.unsafe_wrap(Array,
                                                 pointer(_inverse_transform_matrix),
                                                 (NDIMS, NDIMS,
                                                  ntuple(_ -> nnodes(basis),
                                                         NDIMS)...,
                                                  nelements))

    _jacobian = Vector{RealT}(undef, nnodes(basis)^NDIMS * nelements)
    jacobian = Trixi.unsafe_wrap(Array, pointer(_jacobian),
                                 (ntuple(_ -> nnodes(basis), NDIMS)...,
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

    elements = P4estElementContainerCovariant{NDIMS, RealT, uEltype, NDIMS + 1,
                                              NDIMS + 2,
                                              NDIMS + 3}(node_coordinates,
                                                         jacobian_matrix,
                                                         transform_matrix,
                                                         inverse_transform_matrix,
                                                         jacobian,
                                                         inverse_jacobian,
                                                         surface_flux_values,
                                                         _node_coordinates,
                                                         _jacobian_matrix,
                                                         _transform_matrix,
                                                         _inverse_transform_matrix,
                                                         _jacobian,
                                                         _inverse_jacobian,
                                                         _surface_flux_values)

    init_elements_2d_manifold_in_3d!(elements, mesh, equations, basis)

    return elements
end

# Compute the node positions and metric terms for the covariant form, assuming that the
# domain is a spherical shell. We do not make any assumptions on the mesh topology (can be # any quadrilateral mesh) and make no approximations aside from the initial polynomial
# interpolation to get the mapping from reference LGL nodes to Cartesian coordinates.
function init_elements_2d_manifold_in_3d!(elements, mesh::P4estMesh{2},
                                          ::AbstractCovariantEquations{2},
                                          basis::LobattoLegendreBasis)
    (; node_coordinates, jacobian_matrix, transform_matrix, inverse_transform_matrix,
    jacobian, inverse_jacobian) = elements

    # Place LGL nodes on the mesh according to the polynomial mapping from reference
    # coordinates (xi1, xi2) to Cartesian coordinates (x1, x2, x3) 
    calc_node_coordinates_2d_shell!(node_coordinates, mesh, basis)

    for element in 1:Trixi.ncells(mesh)

        # Differentiate Cartesian components (x1, x2, x3) with respect to reference 
        # coordinates (xi1, xi2) using the polynomial derivative operator.
        Trixi.calc_jacobian_matrix!(jacobian_matrix, element, node_coordinates, basis)

        for j in eachnode(basis), i in eachnode(basis)

            # Convert from Cartesian components (x1, x2, x3) to spherical coordinates
            # (longitude, latitude, radius), with latitude and longitude in radians
            r = sqrt(node_coordinates[1, i, j, element]^2 +
                     node_coordinates[2, i, j, element]^2 +
                     node_coordinates[3, i, j, element]^2)
            lon = atan(node_coordinates[2, i, j, element],
                       node_coordinates[1, i, j, element])
            lat = asin(node_coordinates[3, i, j, element] / r)

            # Precompute trigonometric terms
            sinlon, coslon = sincos(lon)
            sinlat, coslat = sincos(lat)

            # Differentiate (x1, x2, x3) with respect to (longitude, latitude, radius)
            jacobian_cartesian_wrt_spherical = SMatrix{3, 3}(-r * sinlon * coslat,
                                                             r * coslon * coslat, 0,
                                                             -r * coslon * sinlat,
                                                             -r * sinlon * sinlat,
                                                             r * coslat,
                                                             coslon * coslat,
                                                             sinlon * coslat, sinlat)

            # Use the chain rule to compute the Jacobian matrix of 
            # (longitude, latitude, radius) with respect to (xi1, xi2)
            jacobian_spherical_wrt_reference = jacobian_cartesian_wrt_spherical \
                                               SMatrix{3, 2}(jacobian_matrix[:, :, i, j,
                                                                             element])

            # Calculate the transformation matrix A such that G = AᵀA is the matrix of 
            # covariant metric tensor components. This can also be interpreted as the 
            # components of the covariant basis vectors with respect to the spherical 
            # coordinate system. For details, see Appendix B of the following paper:
            #   Nair, R. D., Thomas, S. J., and Loft, R. D. (2005). 
            #   A discontinuous Galerkin transport scheme on the cubed sphere. 
            #   Monthly Weather Review 133, 814-828.
            transform_matrix[1, 1, i, j, element] = r * coslat *
                                                    jacobian_spherical_wrt_reference[1,
                                                                                     1]
            transform_matrix[2, 1, i, j, element] = r *
                                                    jacobian_spherical_wrt_reference[2,
                                                                                     1]
            transform_matrix[1, 2, i, j, element] = r * coslat *
                                                    jacobian_spherical_wrt_reference[1,
                                                                                     2]
            transform_matrix[2, 2, i, j, element] = r *
                                                    jacobian_spherical_wrt_reference[2,
                                                                                     2]

            # Analytically compute the scalar factor (det(GᵀG)) = det(A)
            jacobian[i, j, element] = transform_matrix[1, 1, i, j, element] *
                                      transform_matrix[2, 2, i, j, element] -
                                      transform_matrix[1, 2, i, j, element] *
                                      transform_matrix[2, 1, i, j, element]
            inverse_jacobian[i, j, element] = 1 / jacobian[i, j, element]

            # Analytically invert the transformation matrix
            inverse_transform_matrix[1, 1, i, j, element] = inverse_jacobian[i, j,
                                                                             element] *
                                                            transform_matrix[2, 2, i, j,
                                                                             element]
            inverse_transform_matrix[1, 2, i, j, element] = -inverse_jacobian[i, j,
                                                                              element] *
                                                            transform_matrix[1, 2, i, j,
                                                                             element]
            inverse_transform_matrix[2, 1, i, j, element] = -inverse_jacobian[i, j,
                                                                              element] *
                                                            transform_matrix[2, 1, i, j,
                                                                             element]
            inverse_transform_matrix[2, 2, i, j, element] = inverse_jacobian[i, j,
                                                                             element] *
                                                            transform_matrix[1, 1, i, j,
                                                                             element]
        end
    end

    return nothing
end
end # @muladd
