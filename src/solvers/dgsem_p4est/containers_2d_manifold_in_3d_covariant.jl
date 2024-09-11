@muladd begin
#! format: noindent

# Specialized element container for partial differential equations in covariant form
mutable struct P4estElementContainerCovariant{NDIMS, RealT <: Real, uEltype <: Real,
                                              NDIMSP1, NDIMSP2,
                                              NDIMSP3} <: Trixi.AbstractContainer
    # Physical Cartesian coordinates at each node
    node_coordinates::Array{RealT, NDIMSP2}  # [orientation, node_i, node_j, node_k, element]
    # Transformation matrix from contravariant to zonal/meridional vector components
    transform_matrix::Array{RealT, NDIMSP3}  # [:, :, i, node_i, node_j, node_k, element]
    # Determinant of the transform matrix
    jacobian::Array{RealT, NDIMSP1} # [node_i, node_j, node_k, element]
    inverse_jacobian::Array{RealT, NDIMSP1}  # [node_i, node_j, node_k, element]
    inverse_transform::Array{RealT, NDIMSP3}  # [:, :, i, node_i, node_j, node_k, element]
    # Buffer for calculated surface flux
    surface_flux_values::Array{uEltype, NDIMSP2} # [variable, i, j, direction, element]

    # internal `resize!`able storage
    _node_coordinates::Vector{RealT}
    _transform_matrix::Vector{RealT}
    _jacobian::Vector{RealT}
    _inverse_jacobian::Vector{RealT}
    _inverse_transform::Vector{RealT}
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

@inline function Trixi.get_node_coords(x, ::AbstractCovariantEquations{2}, ::DG,
                                       indices...)
    return SVector(ntuple(@inline(idx->x[idx, indices...]), 3))
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

    _transform_matrix = Vector{RealT}(undef,
                                      NDIMS * NDIMS * nnodes(basis)^NDIMS *
                                      nelements)
    transform_matrix = Trixi.unsafe_wrap(Array, pointer(_transform_matrix),
                                         (NDIMS, NDIMS,
                                          ntuple(_ -> nnodes(basis), NDIMS)...,
                                          nelements))

    _jacobian = Vector{RealT}(undef, nnodes(basis)^NDIMS * nelements)
    jacobian = Trixi.unsafe_wrap(Array, pointer(_jacobian),
                                 (ntuple(_ -> nnodes(basis), NDIMS)...,
                                  nelements))

    _inverse_jacobian = Vector{RealT}(undef, nnodes(basis)^NDIMS * nelements)
    inverse_jacobian = Trixi.unsafe_wrap(Array, pointer(_inverse_jacobian),
                                         (ntuple(_ -> nnodes(basis), NDIMS)...,
                                          nelements))

    _inverse_transform = Vector{RealT}(undef,
                                       NDIMS * NDIMS * nnodes(basis)^NDIMS *
                                       nelements)
    inverse_transform = Trixi.unsafe_wrap(Array,
                                          pointer(_inverse_transform),
                                          (NDIMS, NDIMS,
                                           ntuple(_ -> nnodes(basis),
                                                  NDIMS)...,
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
                                                         transform_matrix,
                                                         jacobian,
                                                         inverse_jacobian,
                                                         inverse_transform,
                                                         surface_flux_values,
                                                         _node_coordinates,
                                                         _transform_matrix,
                                                         _jacobian,
                                                         _inverse_jacobian,
                                                         _inverse_transform,
                                                         _surface_flux_values)

    init_elements_2d_manifold_in_3d!(elements, mesh, equations, basis)

    return elements
end

# Compute the node positions and metric terms for the covariant form, assuming that the
# domain is a spherical shell. We do not make any assumptions on the mesh topology.
function init_elements_2d_manifold_in_3d!(elements, mesh::P4estMesh{2},
                                          equations::AbstractCovariantEquations{2},
                                          basis::LobattoLegendreBasis)
    (; node_coordinates, transform_matrix, inverse_transform,
    jacobian, inverse_jacobian) = elements

    # The tree node coordinates are assumed to be on the spherical shell centred around the # origin. Therefore, we can compute the radius once and use it throughout.
    radius = sqrt(mesh.tree_node_coordinates[1, 1, 1, 1]^2 +
                  mesh.tree_node_coordinates[2, 1, 1, 1]^2 +
                  mesh.tree_node_coordinates[3, 1, 1, 1]^2)

    for element in 1:Trixi.ncells(mesh)

        # For now, we will only use the corner nodes from the P4estMesh, and construct the
        # mapping and its associated metric terms analytically without any polynomial 
        # approximation following the approach described in by Guba et al. (see 
        # https://doi.org/10.5194/gmd-7-2803-2014, Appendix A). If the option 
        # "element_local_mapping = true" is used when constructing the mesh using
        # P4estMeshCubedSphere2D and the polynomial degree of the mesh is the same as that 
        # of the solver, the node positions are equal to those obtained using the standard
        # calc_node_coordinates! function.

        # Extract the corner vertex positions from the P4estMesh
        nnodes = length(mesh.nodes)
        v1 = SVector{3}(view(mesh.tree_node_coordinates, :, 1, 1, element))
        v2 = SVector{3}(view(mesh.tree_node_coordinates, :, nnodes, 1, element))
        v3 = SVector{3}(view(mesh.tree_node_coordinates, :, nnodes, nnodes, element))
        v4 = SVector{3}(view(mesh.tree_node_coordinates, :, 1, nnodes, element))

        for j in eachnode(basis), i in eachnode(basis)

            # get reference node coordinates
            xi1, xi2 = basis.nodes[i], basis.nodes[j]

            # Construct a bilinear mapping based on the four corner vertices
            x_bilinear = 0.25f0 *
                         ((1 - xi1) * (1 - xi2) * v1 + (1 + xi1) * (1 - xi2) * v2 +
                          (1 + xi1) * (1 + xi2) * v3 + (1 - xi1) * (1 + xi2) * v4)

            # Project the mapped local coordinates onto the sphere 
            scaling_factor = radius / norm(x_bilinear)
            x = scaling_factor * x_bilinear

            # Convert to longitude and latitude
            lon, lat = atan(x[2], x[1]), asin(x[3] / radius)

            # Compute trigonometric terms needed for analytical metrics
            sinlon, coslon = sincos(lon)
            sinlat, coslat = sincos(lat)
            a11 = sinlon * sinlon * coslat * coslat + sinlat * sinlat
            a12 = -sinlon * coslon * coslat * coslat
            a13 = -coslon * sinlat * coslat
            a21 = a12
            a22 = coslon * coslon * coslat * coslat + sinlat * sinlat
            a23 = -sinlon * sinlat * coslat
            a31 = -coslon * sinlat
            a32 = -sinlon * sinlat
            a33 = coslat

            # Analytically compute the transformation matrix A, such that G = AᵀA is the 
            # covariant metric tensor and a_i = A[1,i] * e_lon + A[2,i] * e_lat denotes 
            # the covariant tangent basis, where e_lon and e_lat are the unit basis vectors
            # in the longitudinal and latitudinal directions, respectively.
            A = 0.25f0 * scaling_factor * SMatrix{2, 3}(-sinlon, 0, coslon, 0, 0, 1) *
                SMatrix{3, 3}(a11, a21, a31, a12, a22, a32, a13, a23, a33) *
                SMatrix{3, 4}(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3],
                              v3[1], v3[2], v3[3], v4[1], v4[2], v4[3]) *
                SMatrix{4, 2}(-1 + xi2, 1 - xi2, 1 + xi2, -1 - xi2,
                              -1 + xi1, -1 - xi1, 1 + xi1, 1 - xi1)

            # Set variables in the cache
            node_coordinates[1, i, j, element] = x[1]
            node_coordinates[2, i, j, element] = x[2]
            node_coordinates[3, i, j, element] = x[3]

            transform_matrix[1, 1, i, j, element] = A[1, 1]
            transform_matrix[2, 1, i, j, element] = A[2, 1]
            transform_matrix[1, 2, i, j, element] = A[1, 2]
            transform_matrix[2, 2, i, j, element] = A[2, 2]

            jacobian[i, j, element] = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
            inverse_jacobian[i, j, element] = 1 / jacobian[i, j, element]

            inverse_transform[1, 1, i, j, element] = inverse_jacobian[i, j, element] *
                                                     A[2, 2]
            inverse_transform[1, 2, i, j, element] = -inverse_jacobian[i, j, element] *
                                                     A[1, 2]
            inverse_transform[2, 1, i, j, element] = -inverse_jacobian[i, j, element] *
                                                     A[2, 1]
            inverse_transform[2, 2, i, j, element] = inverse_jacobian[i, j, element] *
                                                     A[1, 1]
        end
    end

    return nothing
end
end # @muladd
