@muladd begin
#! format: noindent

# Specialized element container for partial differential equations in covariant form
mutable struct P4estElementContainerCovariant{NDIMS, RealT <: Real, uEltype <: Real,
                                              NDIMSP1, NDIMSP2,
                                              NDIMSP3, NDIMSP4} <:
               Trixi.AbstractContainer
    # Physical Cartesian coordinates at each node
    # [orientation, node_i, node_j, node_k, element]
    node_coordinates::Array{RealT, NDIMSP2}
    # Covariant basis for the tangent space, expanded with respect to the spherical basis
    # [component index, basis vector index, i, node_i, node_j, node_k, element]
    covariant_basis::Array{RealT, NDIMSP3}
    # Covariant components of the metric tensor G_{lm}
    covariant_metric::Array{RealT, NDIMSP3}  # [l, m, node_i, node_j, node_k, element]
    # Contravariant components of the metric tensor G^{lm}
    contravariant_metric::Array{RealT, NDIMSP3}  # [l, m, node_i, node_j, node_k, element]
    # Christoffel symbols of the second kind Γ_{lm}^k
    christoffel_symbols::Array{RealT, NDIMSP4}  # [k, l, m, node_i, node_j, node_k, element]
    inverse_jacobian::Array{RealT, NDIMSP1}  # [node_i, node_j, node_k, element]
    surface_flux_values::Array{uEltype, NDIMSP2} # [variable, i, j, direction, element]

    # internal `resize!`able storage
    _node_coordinates::Vector{RealT}
    _covariant_basis::Vector{RealT}
    _covariant_metric::Vector{RealT}
    _contravariant_metric::Vector{RealT}
    _christoffel_symbols::Vector{RealT}
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

    _christoffel_symbols = Vector{RealT}(undef,
                                         NDIMS * NDIMS * NDIMS *
                                         nnodes(basis)^NDIMS * nelements)

    christoffel_symbols = Trixi.unsafe_wrap(Array, pointer(_christoffel_symbols),
                                            (NDIMS, NDIMS, NDIMS,
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
                                              NDIMS + 3,
                                              NDIMS + 4}(node_coordinates,
                                                         covariant_basis,
                                                         covariant_metric,
                                                         contravariant_metric,
                                                         christoffel_symbols,
                                                         inverse_jacobian,
                                                         surface_flux_values,
                                                         _node_coordinates,
                                                         _covariant_basis,
                                                         _covariant_metric,
                                                         _contravariant_metric,
                                                         _christoffel_symbols,
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
    (; node_coordinates, covariant_basis, covariant_metric, contravariant_metric,
    christoffel_symbols, inverse_jacobian) = elements

    # The tree node coordinates are assumed to be on the spherical shell centred around the
    # origin. Therefore, we can compute the radius once and use it throughout.
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

            # Get reference node coordinates
            xi1, xi2 = basis.nodes[i], basis.nodes[j]

            # Get Cartesian coordinates and covariant basis using formulas from Guba et al.
            x = local_mapping(xi1, xi2, v1, v2, v3, v4, radius)
            A_cov = local_covariant_basis(xi1, xi2, v1, v2, v3, v4, radius)

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

            inverse_jacobian[i, j, element] = 1 / sqrt(G_cov[1, 1] * G_cov[2, 2] -
                                                   G_cov[1, 2] * G_cov[2, 1])
        end

        calc_christoffel_symbols!(christoffel_symbols, covariant_basis,
                                  contravariant_metric,
                                  mesh, basis)
    end

    return nothing
end

# Approximate the Christoffel symbols by differentiating the (orthonormal) 
# spherical components of the covariant basis vectors with respect to the reference 
# coordinates (using the discrete derivative operator) and taking the dot product with the 
# contravariant basis to obtain Γ_{l,m}^k = Γ_{m,l}^k ≈ (aᵏ)ₙ⋅∂Iᴺ(aₘ)ₙ/∂ξˡ 
function calc_christoffel_symbols!(christoffel_symbols, covariant_basis,
                                   contravariant_metric, mesh::P4estMesh{2},
                                   basis::LobattoLegendreBasis)
    (; derivative_matrix) = basis

    for element in 1:Trixi.ncells(mesh)
        for j in eachnode(basis), i in eachnode(basis)
            for k in 1:2
                # Raise indices of the covariant basis to obtain contravariant basis
                a_con_1 = contravariant_metric[k, 1, i, j, element] *   # (aᵏ)₁
                          covariant_basis[1, 1, i, j, element] +
                          contravariant_metric[k, 2, i, j, element] *
                          covariant_basis[1, 2, i, j, element]
                a_con_2 = contravariant_metric[k, 1, i, j, element] *  # (aᵏ)₂
                          covariant_basis[2, 1, i, j, element] +
                          contravariant_metric[k, 2, i, j, element] *
                          covariant_basis[2, 2, i, j, element]
                for m in 1:2
                    # Discretely differentiate covariant basis components with respect to ξ¹
                    da1dxi1 = zero(eltype(christoffel_symbols))  # ∂Iᴺ(aₘ)₁/∂ξ¹
                    da2dxi1 = zero(eltype(christoffel_symbols))  # ∂Iᴺ(aₘ)₂/∂ξ¹
                    for ii in eachnode(basis)
                        da1dxi1 = da1dxi1 +
                                  derivative_matrix[ii, j] *
                                  covariant_basis[1, m, ii, j, element]
                        da2dxi1 = da2dxi1 +
                                  derivative_matrix[ii, j] *
                                  covariant_basis[2, m, ii, j, element]
                    end
                    # Discretely differentiate covariant basis components with respect to ξ²
                    da1dxi2 = zero(eltype(christoffel_symbols))  # ∂Iᴺ(aₘ)₁/∂ξ²
                    da2dxi2 = zero(eltype(christoffel_symbols))  # ∂Iᴺ(aₘ)₂/∂ξ²
                    for jj in eachnode(basis)
                        da1dxi2 = da1dxi2 +
                                  derivative_matrix[i, jj] *
                                  covariant_basis[1, m, i, jj, element]
                        da2dxi2 = da2dxi2 +
                                  derivative_matrix[i, jj] *
                                  covariant_basis[2, m, i, jj, element]
                    end
                    # Take the dot product (assuming an orthonormal basis)
                    # Γ_{m,1}^k = (aᵏ)₁⋅∂Iᴺ(aₘ)₁/∂ξ¹ + (aᵏ)₂⋅∂Iᴺ(aₘ)₂/∂ξ¹
                    christoffel_symbols[k, m, 1, i, j, element] = a_con_1 * da1dxi1 +
                                                                  a_con_2 * da2dxi1
                    # Γ_{m,2}^k = (aᵏ)₁⋅∂Iᴺ(aₘ)₁/∂ξ² + (aᵏ)₂⋅∂Iᴺ(aₘ)₂/∂ξ²
                    christoffel_symbols[k, m, 2, i, j, element] = a_con_1 * da1dxi2 +
                                                                  a_con_2 * da2dxi2
                end
            end
        end
    end
end
end # @muladd
