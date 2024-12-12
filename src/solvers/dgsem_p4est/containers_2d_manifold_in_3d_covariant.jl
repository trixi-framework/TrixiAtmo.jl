@muladd begin
#! format: noindent

# Container for storing values of auxiliary variables at volume/surface quadrature nodes
struct P4estAuxiliaryNodeVariableContainer{NDIMS, uEltype <: Real, NDIMSP2}
    aux_node_vars::Array{uEltype, NDIMSP2} # [variable, i, j, k, element]
    aux_surface_node_vars::Array{uEltype, NDIMSP2} # [variable, i, j, k, element]

    # internal `resize!`able storage
    _aux_node_vars::Vector{uEltype}
    _aux_surface_node_vars::Vector{uEltype}
end

# Return the auxiliary variables at a given volume node index
@inline function get_node_aux_vars(aux_node_vars, equations, ::DG, indices...)
    return SVector(ntuple(@inline(v->aux_node_vars[v, indices...]),
                          Val(n_aux_node_vars(equations))))
end

# Return the auxiliary variables at a given surface node index
@inline function get_surface_node_aux_vars(aux_surface_node_vars, equations, ::DG,
                                           indices...)
    aux_vars_ll = SVector(ntuple(@inline(v->aux_surface_node_vars[1, v, indices...]),
                                 Val(n_aux_node_vars(equations))))
    aux_vars_rr = SVector(ntuple(@inline(v->aux_surface_node_vars[2, v, indices...]),
                                 Val(n_aux_node_vars(equations))))
    return aux_vars_ll, aux_vars_rr
end

# Create auxiliary node variable container and initialize auxiliary variables
function init_auxiliary_node_variables(mesh::Union{P4estMesh, T8codeMesh},
                                       equations, dg, elements, interfaces)
    nelements = Trixi.ncells(mesh)
    ninterfaces = Trixi.count_required_surfaces(mesh).interfaces
    NDIMS = ndims(elements)
    uEltype = eltype(elements)

    _aux_node_vars = Vector{uEltype}(undef,
                                     n_aux_node_vars(equations) *
                                     nnodes(dg)^NDIMS * nelements)
    aux_node_vars = Trixi.unsafe_wrap(Array, pointer(_aux_node_vars),
                                      (n_aux_node_vars(equations),
                                       ntuple(_ -> nnodes(dg), NDIMS)...,
                                       nelements))
    _aux_surface_node_vars = Vector{uEltype}(undef,
                                             2 * n_aux_node_vars(equations) *
                                             nnodes(dg)^(NDIMS - 1) *
                                             ninterfaces)
    aux_surface_node_vars = Trixi.unsafe_wrap(Array,
                                              pointer(_aux_surface_node_vars),
                                              (2, n_aux_node_vars(equations),
                                               ntuple(_ -> nnodes(dg),
                                                      NDIMS - 1)...,
                                               ninterfaces))

    auxiliary_variables = P4estAuxiliaryNodeVariableContainer{NDIMS,
                                                              uEltype,
                                                              NDIMS + 2}(aux_node_vars,
                                                                         aux_surface_node_vars,
                                                                         _aux_node_vars,
                                                                         _aux_surface_node_vars)

    init_auxiliary_node_variables!(auxiliary_variables, mesh, equations, dg, elements)
    init_auxiliary_surface_node_variables!(auxiliary_variables, mesh, equations, dg,
                                           interfaces)
    return auxiliary_variables
end

# By default, the auxiliary surface node variables are just the volume node variables 
# evaluated at the surface nodes, similarly to prolong2interfaces. This could, however, be 
# specialized for other applications.
function init_auxiliary_surface_node_variables!(auxiliary_variables::P4estAuxiliaryNodeVariableContainer{2},
                                                mesh::Union{P4estMesh{2},
                                                            T8codeMesh{2}}, equations,
                                                dg::DG, interfaces)
    (; aux_node_vars, aux_surface_node_vars) = auxiliary_variables
    index_range = eachnode(dg)

    Trixi.@threaded for interface in axes(interfaces.node_indices, 2)
        # Copy solution data from the primary element using "delayed indexing" with
        # a start value and a step size to get the correct face and orientation.
        # Note that in the current implementation, the interface will be
        # "aligned at the primary element", i.e., the index of the primary side
        # will always run forwards.
        primary_element = interfaces.neighbor_ids[1, interface]
        primary_indices = interfaces.node_indices[1, interface]

        i_primary_start, i_primary_step = Trixi.index_to_start_step_2d(primary_indices[1],
                                                                       index_range)
        j_primary_start, j_primary_step = Trixi.index_to_start_step_2d(primary_indices[2],
                                                                       index_range)

        i_primary = i_primary_start
        j_primary = j_primary_start
        for i in eachnode(dg)
            for v in axes(aux_node_vars, 1)
                aux_surface_node_vars[1, v, i, interface] = aux_node_vars[v,
                                                                          i_primary,
                                                                          j_primary,
                                                                          primary_element]
            end
            i_primary += i_primary_step
            j_primary += j_primary_step
        end

        # Copy solution data from the secondary element using "delayed indexing" with
        # a start value and a step size to get the correct face and orientation.
        secondary_element = interfaces.neighbor_ids[2, interface]
        secondary_indices = interfaces.node_indices[2, interface]

        i_secondary_start, i_secondary_step = Trixi.index_to_start_step_2d(secondary_indices[1],
                                                                           index_range)
        j_secondary_start, j_secondary_step = Trixi.index_to_start_step_2d(secondary_indices[2],
                                                                           index_range)

        i_secondary = i_secondary_start
        j_secondary = j_secondary_start
        for i in eachnode(dg)
            for v in axes(aux_node_vars, 1)
                aux_surface_node_vars[2, v, i, interface] = aux_node_vars[v,
                                                                          i_secondary,
                                                                          j_secondary,
                                                                          secondary_element]
            end
            i_secondary += i_secondary_step
            j_secondary += j_secondary_step
        end
    end

    return nothing
end

# Get Cartesian node positions for the covariant form, dispatching on the dimension of the 
# manifold as well as the ambient dimension
@inline function Trixi.get_node_coords(x,
                                       ::AbstractCovariantEquations{NDIMS,
                                                                    NDIMS_AMBIENT},
                                       ::DG,
                                       indices...) where {NDIMS, NDIMS_AMBIENT}
    return SVector(ntuple(@inline(idx->x[idx, indices...]), NDIMS_AMBIENT))
end

# Compute the auxiliary metric terms for the covariant form, assuming that the
# domain is a spherical shell. We do not make any assumptions on the mesh topology, but we 
# require that the elements are constructed using the element-local mapping from the 
# following paper:
#
# O. Guba, M. A. Taylor, P. A. Ullrich, J. R. Overfelt, and M. N. Levy (2014). The spectral
# element method (SEM) on variable-resolution grids: evaluating grid sensitivity and 
# resolution-aware numerical viscosity. Geoscientific Model Development 7(6) 2803–2816.
# DOI: 10.5194/gmd-7-2803-2014
# 
# When creating a cubed sphere mesh using P4estMeshCubedSphere2D, this is enabled by 
# passing the keyword argument "element_local_mapping = true" to the constructor and 
# ensuring that the polynomial degree of the mesh is equal to that of the solver.
#
# Otherwise, the mesh's tree node coordinates will be interpolated to the solver's 
# physical node coordinates, and this would introduce a polynomial approximation of the 
# geometry, making the analytical metric terms computed here no longer correct.
function init_auxiliary_node_variables!(auxiliary_variables, mesh::P4estMesh{2, 3},
                                        equations::AbstractCovariantEquations{2, 3}, dg,
                                        elements)
    (; tree_node_coordinates) = mesh
    (; aux_node_vars) = auxiliary_variables

    NDIMS = 2
    NDIMS_AMBIENT = 3

    # Check that the degree of the mesh matches that of the solver
    @assert length(mesh.nodes) == nnodes(dg)

    # The tree node coordinates are assumed to be on the spherical shell centred around the 
    # origin. Therefore, we can compute the radius once and use it throughout.
    radius = norm(Trixi.get_node_coords(tree_node_coordinates, equations, dg, 1, 1, 1))

    Trixi.@threaded for element in 1:Trixi.ncells(mesh)
        # Extract the corner vertex positions from the P4estMesh
        v1 = Trixi.get_node_coords(tree_node_coordinates, equations, dg,
                                   1, 1, element)
        v2 = Trixi.get_node_coords(tree_node_coordinates, equations, dg,
                                   nnodes(dg), 1, element)
        v3 = Trixi.get_node_coords(tree_node_coordinates, equations, dg,
                                   nnodes(dg), nnodes(dg), element)
        v4 = Trixi.get_node_coords(tree_node_coordinates, equations, dg,
                                   1, nnodes(dg), element)

        # Compute the auxiliary metric information at each node
        for j in eachnode(dg), i in eachnode(dg)
           

            # Covariant basis in the desired global coordinate system as columns of a matrix
            basis_covariant = calc_basis_covariant(v1, v2, v3, v4,
                                    dg.basis.nodes[i], dg.basis.nodes[j],
                                    radius,
                                    equations.global_coordinate_system)
            start_ind, end_ind = 1, length(basis_covariant)
            aux_node_vars[start_ind:end_ind, i, j, element] = SVector(basis_covariant)
            start_ind = end_ind + 1

            # Metric tensor 
            metric_covariant = basis_covariant' * basis_covariant
            metric_contravariant = inv(basis_covariant' * basis_covariant)

            # Contravariant basis vectors as rows of a matrix
            basis_contravariant = metric_contravariant * basis_covariant'
            end_ind = start_ind+length(basis_contravariant)-1
            aux_node_vars[start_ind:end_ind, i, j, element] = SVector(basis_contravariant)
            start_ind = end_ind + 1

            # Area element
            area_element = sqrt(det(metric_covariant))
            aux_node_vars[start_ind, i, j, element] = area_element
            start_ind = start_ind + 1

            # Covariant metric tensor components
            end_ind = start_ind + 2
            aux_node_vars[start_ind:end_ind, i, j, element] = SVector(metric_covariant[1,1], metric_covariant[1,2], metric_covariant[2,2])
            start_ind = end_ind + 1

            # Contravariant metric tensor components
            end_ind = start_ind + 2
            aux_node_vars[start_ind:end_ind, i, j, element] = SVector(metric_contravariant[1,1], metric_contravariant[1,2], metric_contravariant[2,2])
        end
    end
    return nothing
end

# Analytically compute the transformation matrix A, such that G = AᵀA is the 
# covariant metric tensor and a_i = A[1,i] * e_x + A[2,i] * e_y + A[3,i] * e_z denotes 
# the covariant tangent basis, where e_x, e_y, and e_z are the Cartesian unit basis vectors.
@inline function calc_basis_covariant(v1, v2, v3, v4, xi1, xi2, radius,
                                      ::GlobalCartesianCoordinates)

    # Construct a bilinear mapping based on the four corner vertices
    xe = 0.25f0 * ((1 - xi1) * (1 - xi2) * v1 + (1 + xi1) * (1 - xi2) * v2 +
          (1 + xi1) * (1 + xi2) * v3 + (1 - xi1) * (1 + xi2) * v4)

    # Derivatives of bilinear map with respect to reference coordinates xi1, xi2
    dxedxi1 = 0.25f0 *
              (-(1 - xi2) * v1 + (1 - xi2) * v2 + (1 + xi2) * v3 - (1 + xi2) * v4)
    dxedxi2 = 0.25f0 *
              (-(1 - xi1) * v1 - (1 + xi1) * v2 + (1 + xi1) * v3 + (1 - xi1) * v4)

    # Use product/quotient rule on the projection
    norm_xe = norm(xe)
    dxdxi1 = radius / norm_xe * (dxedxi1 - dot(xe, dxedxi1) / norm_xe^2 * xe)
    dxdxi2 = radius / norm_xe * (dxedxi2 - dot(xe, dxedxi2) / norm_xe^2 * xe)

    return SMatrix{3, 2}(dxdxi1[1], dxdxi1[2], dxdxi1[3],
                         dxdxi2[1], dxdxi2[2], dxdxi2[3])
end

# Analytically compute the transformation matrix A, such that G = AᵀA is the 
# covariant metric tensor and a_i = A[1,i] * e_lon + A[2,i] * e_lat denotes 
# the covariant tangent basis, where e_lon and e_lat are the unit basis vectors
# in the longitudinal and latitudinal directions, respectively. This formula is 
# taken from Guba et al. (2014).
@inline function calc_basis_covariant(v1, v2, v3, v4, xi1, xi2, radius,
                                      ::GlobalSphericalCoordinates)
    # Construct a bilinear mapping based on the four corner vertices
    xe = 0.25f0 * ((1 - xi1) * (1 - xi2) * v1 + (1 + xi1) * (1 - xi2) * v2 +
          (1 + xi1) * (1 + xi2) * v3 + (1 - xi1) * (1 + xi2) * v4)

    # Project the mapped local coordinates onto the sphere using a simple scaling
    scaling_factor = radius / norm(xe)
    x = scaling_factor * xe

    # Convert Cartesian coordinates to longitude and latitude
    lon, lat = atan(x[2], x[1]), asin(x[3] / radius)

    # Compute trigonometric terms needed for analytical metrics
    sinlon, coslon = sincos(lon)
    sinlat, coslat = sincos(lat)
    a11 = sinlon * sinlon * coslat * coslat + sinlat * sinlat
    a12 = a21 = -sinlon * coslon * coslat * coslat
    a13 = -coslon * sinlat * coslat
    a22 = coslon * coslon * coslat * coslat + sinlat * sinlat
    a23 = -sinlon * sinlat * coslat
    a31 = -coslon * sinlat
    a32 = -sinlon * sinlat
    a33 = coslat

    # Compute the matrix A containing spherical components of the covariant basis
    A = 0.25f0 * scaling_factor *
        SMatrix{2, 3}(-sinlon, 0, coslon, 0, 0, 1) *
        SMatrix{3, 3}(a11, a21, a31, a12, a22, a32, a13, a23, a33) *
        SMatrix{3, 4}(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3],
                      v3[1], v3[2], v3[3], v4[1], v4[2], v4[3]) *
        SMatrix{4, 2}(-1 + xi2, 1 - xi2, 1 + xi2, -1 - xi2,
                      -1 + xi1, -1 - xi1, 1 + xi1, 1 - xi1)

    # Make zero component in the radial direction so the matrix has the right dimensions
    return SMatrix{3, 2}(A[1, 1], A[2, 1], 0.0f0, A[1, 2], A[2, 2], 0.0f0)
end
end # muladd
