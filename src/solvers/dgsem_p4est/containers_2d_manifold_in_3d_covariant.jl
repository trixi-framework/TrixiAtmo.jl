@muladd begin
#! format: noindent

# Specialized element container for partial differential equations in covariant form
mutable struct P4estElementContainerCovariant{NDIMS, RealT <: Real, uEltype <: Real,
                                              NDIMSP1, NDIMSP2,
                                              NDIMSP3} <: Trixi.AbstractContainer
    # Physical Cartesian coordinates at each node
    # [orientation, node_i, node_j, node_k, element]
    node_coordinates::Array{RealT, NDIMSP2}
    # Covariant basis for the tangent space, expanded with respect to the spherical basis
    covariant_basis::Array{RealT, NDIMSP3}  # [:, :, i, node_i, node_j, node_k, element]
    covariant_metric::Array{RealT, NDIMSP3}
    contravariant_metric::Array{RealT, NDIMSP3}
    inverse_jacobian::Array{RealT, NDIMSP1}  # [node_i, node_j, node_k, element]
    surface_flux_values::Array{uEltype, NDIMSP2} # [variable, i, j, direction, element]

    # internal `resize!`able storage
    _node_coordinates::Vector{RealT}
    _covariant_basis::Vector{RealT}
    _covariant_metric::Vector{RealT}
    _contravariant_metric::Vector{RealT}
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

# Get Cartesian node positions, dispatching on the dimension of the equation system
@inline function Trixi.get_node_coords(x,
                                       ::AbstractCovariantEquations{NDIMS,
                                                                    NDIMS_AMBIENT},
                                       ::DG,
                                       indices...) where {NDIMS, NDIMS_AMBIENT}
    return SVector(ntuple(@inline(idx->x[idx, indices...]), NDIMS_AMBIENT))
end

# Volume element J = sqrt(det(G)), where G is the matrix of covariant metric components
@inline function volume_element(elements, i, j, element)
    return 1 / elements.inverse_jacobian[i, j, element]
end

# Convert contravariant velocity/momentum components to zonal and meridional components
@inline function contravariant2spherical(v_con_1, v_con_2, elements, i, j, element)
    return elements.covariant_basis[1, 1, i, j, element] * v_con_1 +
           elements.covariant_basis[1, 2, i, j, element] * v_con_2,
           elements.covariant_basis[2, 1, i, j, element] * v_con_1 +
           elements.covariant_basis[2, 2, i, j, element] * v_con_2
end

# Convert zonal and meridional velocity/momentum components to contravariant components
@inline function spherical2contravariant(v_lon, v_lat, elements, i, j, element)
    Jv_con_1 = elements.covariant_basis[2, 2, i, j, element] * v_lon -
               elements.covariant_basis[1, 2, i, j, element] * v_lat
    Jv_con_2 = -elements.covariant_basis[2, 1, i, j, element] * v_lon +
               elements.covariant_basis[1, 1, i, j, element] * v_lat
    return Jv_con_1 * elements.inverse_jacobian[i, j, element],
           Jv_con_2 * elements.inverse_jacobian[i, j, element]
end

# Create element container and initialize element data for a mesh of dimension NDIMS in 
# an ambient space of dimension NDIMS_AMBIENT, with the equations expressed in covariant 
# form
function Trixi.init_elements(mesh::P4estMesh{NDIMS, NDIMS_AMBIENT, RealT},
                             equations::AbstractCovariantEquations{NDIMS,
                                                                   NDIMS_AMBIENT},
                             basis,
                             ::Type{uEltype}) where {NDIMS,
                                                     NDIMS_AMBIENT,
                                                     RealT <: Real,
                                                     uEltype <: Real}
    nelements = Trixi.ncells(mesh)

    _node_coordinates = Vector{RealT}(undef,
                                      NDIMS_AMBIENT * nnodes(basis)^NDIMS * nelements)
    node_coordinates = Trixi.unsafe_wrap(Array, pointer(_node_coordinates),
                                         (NDIMS_AMBIENT,
                                          ntuple(_ -> nnodes(basis), NDIMS)...,
                                          nelements))

    _covariant_basis = Vector{RealT}(undef,
                                     NDIMS * NDIMS * nnodes(basis)^NDIMS *
                                     nelements)
    covariant_basis = Trixi.unsafe_wrap(Array, pointer(_covariant_basis),
                                        (NDIMS, NDIMS,
                                         ntuple(_ -> nnodes(basis), NDIMS)...,
                                         nelements))

    _covariant_metric = Vector{RealT}(undef,
                                      NDIMS * NDIMS * nnodes(basis)^NDIMS *
                                      nelements)
    covariant_metric = Trixi.unsafe_wrap(Array, pointer(_covariant_metric),
                                         (NDIMS, NDIMS,
                                          ntuple(_ -> nnodes(basis), NDIMS)...,
                                          nelements))

    _contravariant_metric = Vector{RealT}(undef,
                                          NDIMS * NDIMS * nnodes(basis)^NDIMS *
                                          nelements)

    contravariant_metric = Trixi.unsafe_wrap(Array, pointer(_contravariant_metric),
                                             (NDIMS, NDIMS,
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

    elements = P4estElementContainerCovariant{NDIMS, RealT, uEltype, NDIMS + 1,
                                              NDIMS + 2,
                                              NDIMS + 3}(node_coordinates,
                                                         covariant_basis,
                                                         covariant_metric,
                                                         contravariant_metric,
                                                         inverse_jacobian,
                                                         surface_flux_values,
                                                         _node_coordinates,
                                                         _covariant_basis,
                                                         _covariant_metric,
                                                         _contravariant_metric,
                                                         _inverse_jacobian,
                                                         _surface_flux_values)

    Trixi.init_elements!(elements, mesh, equations, basis)

    return elements
end

# Compute the node positions and metric terms for the covariant form, assuming that the
# domain is a spherical shell. We do not make any assumptions on the mesh topology.
function Trixi.init_elements!(elements, mesh::P4estMesh{2, 3},
                              equations::AbstractCovariantEquations{2, 3},
                              basis::LobattoLegendreBasis)
    (; node_coordinates, covariant_basis,
    covariant_metric, contravariant_metric, inverse_jacobian) = elements

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
            A_cov = 0.25f0 * scaling_factor *
                    SMatrix{2, 3}(-sinlon, 0, coslon, 0, 0, 1) *
                    SMatrix{3, 3}(a11, a21, a31, a12, a22, a32, a13, a23, a33) *
                    SMatrix{3, 4}(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3],
                                  v3[1], v3[2], v3[3], v4[1], v4[2], v4[3]) *
                    SMatrix{4, 2}(-1 + xi2, 1 - xi2, 1 + xi2, -1 - xi2,
                                  -1 + xi1, -1 - xi1, 1 + xi1, 1 - xi1)

            # Covariant and contravariant metric tensor components
            G_cov = A_cov' * A_cov
            G_con = inv(G_cov)

            # Set variables in the cache
            node_coordinates[1, i, j, element] = x[1]
            node_coordinates[2, i, j, element] = x[2]
            node_coordinates[3, i, j, element] = x[3]

            for l in 1:2, m in 1:2
                covariant_basis[l, m, i, j, element] = A_cov[l, m]
                covariant_metric[l, m, i, j, element] = G_cov[l, m]
                contravariant_metric[l, m, i, j, element] = G_con[l, m]
            end

            inverse_jacobian[i, j, element] = 1 /
                                              sqrt(G_cov[1, 1] * G_cov[2, 2] -
                                                   G_cov[1, 2] * G_cov[2, 1])
        end
    end

    return nothing
end
end # @muladd
