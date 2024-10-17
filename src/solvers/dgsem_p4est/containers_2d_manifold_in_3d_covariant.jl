@muladd begin
#! format: noindent

@doc raw"""
    P4estElementContainerCovariant{NDIMS, RealT <: Real, uEltype <: Real,
                                   NDIMSP1, NDIMSP2,
                                   NDIMSP3, NDIMSP4} <: Trixi.AbstractContainer

Specialized element container for equations in covariant form on a manifold of dimension
`NDIMS`, created when `Trixi.init_elements` is dispatched on 
`AbstractCovariantEquations{NDIMS, NDIMS_AMBIENT}` and `P4estMesh{NDIMS, NDIMS_AMBIENT`.
Contains the following geometric information:
- `node_coordinates::Array{RealT, NDIMSP2}:` Cartesian components of the node positions
  within the ambient space of dimension `NDIMS_AMBIENT`, where the first index is the 
  Cartesian component index, the next `NDIMS` indices are tensor-product node indices ($i$, 
  $j$, for `NDIMS = 2` or $i$, $j$, $k$, for `NDIMS = 3`), and the last is the element
  index.
- `covariant_basis::Array{RealT, NDIMSP3}`: Components of the covariant tangent basis
  vectors  $\vec{a}_m = \partial \vec{X} / \partial \xi^m$ expanded in terms of a global  
  tangent basis (i.e. zonal and meridional components in the case of a spherical shell), 
  where $\vec{X}$ denotes the mapping from the reference element to the physical element.
  he first index is the global component index $l$, the second index is the local
  component index $m$, the next `NDIMS` indices are tensor-product node indices ($i$, $j$, 
  for `NDIMS = 2` or $i$, $j$, $k$, for `NDIMS = 3`), and the last is the element index.
- `inverse_jacobian::Array{RealT, NDIMSP1}`: Nodal values of $1/J$, where 
  $J = \sqrt{\mathrm{det}(\bm{G})}$, and $\bm{G}$ is the matrix of covariant metric tensor 
  components $G_{lm} = \vec{a}_l \cdot \vec{a}_m$. The first `NDIMS` indices are 
  tensor-product node indices ($i$, $j$, for `NDIMS = 2` or $i$, $j$, $k$, for 
  `NDIMS = 3`), and the last is the element index.

!!! warning
    The covariant solver currently only supports conforming meshes.
"""
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
    christoffel_symbols::Array{RealT, NDIMSP4}  # [l, m, k, node_i, node_j, node_k, element]
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

    Trixi.@threaded for element in 1:Trixi.ncells(mesh)

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

        calc_christoffel_symbols!(christoffel_symbols, covariant_metric,
                                  contravariant_metric, mesh, basis)
    end

    return nothing
end

# Calculate Christoffel symbols (this is done approximately using the collocation
# derivative)
function calc_christoffel_symbols!(christoffel_symbols, covariant_metric,
                                   contravariant_metric, mesh::P4estMesh{2},
                                   basis::LobattoLegendreBasis)
    (; derivative_matrix) = basis

    Trixi.@threaded for element in 1:Trixi.ncells(mesh)
        for j in eachnode(basis), i in eachnode(basis)

            # Differentiate covariant metric components with respect to ξ¹
            dG11dxi1 = zero(eltype(christoffel_symbols))
            dG12dxi1 = zero(eltype(christoffel_symbols))
            dG22dxi1 = zero(eltype(christoffel_symbols))
            for ii in eachnode(basis)
                dG11dxi1 = dG11dxi1 +
                           derivative_matrix[i, ii] *
                           covariant_metric[1, 1, ii, j, element]
                dG12dxi1 = dG12dxi1 +
                           derivative_matrix[i, ii] *
                           covariant_metric[1, 2, ii, j, element]
                dG22dxi1 = dG22dxi1 +
                           derivative_matrix[i, ii] *
                           covariant_metric[2, 2, ii, j, element]
            end

            # Differentiate covariant metric components with respect to ξ²
            dG11dxi2 = zero(eltype(christoffel_symbols))
            dG12dxi2 = zero(eltype(christoffel_symbols))
            dG22dxi2 = zero(eltype(christoffel_symbols))
            for jj in eachnode(basis)
                dG11dxi2 = dG11dxi2 +
                           derivative_matrix[j, jj] *
                           covariant_metric[1, 1, i, jj, element]
                dG12dxi2 = dG12dxi2 +
                           derivative_matrix[j, jj] *
                           covariant_metric[1, 2, i, jj, element]
                dG22dxi2 = dG22dxi2 +
                           derivative_matrix[j, jj] *
                           covariant_metric[2, 2, i, jj, element]
            end

            # Compute Christoffel symbols of the first kind
            Gamma_1 = SMatrix{2, 2}(0.5f0 * dG11dxi1, 0.5f0 * dG11dxi2,
                                    0.5f0 * dG11dxi2, dG12dxi2 - 0.5f0 * dG22dxi1)
            Gamma_2 = SMatrix{2, 2}(dG12dxi1 - 0.5f0 * dG11dxi2, 0.5f0 * dG22dxi1,
                                    0.5f0 * dG22dxi1, 0.5f0 * dG22dxi2)

            # Raise indices to get Christoffel symbols of the second kind
            G_con_11 = contravariant_metric[1, 1, i, j, element]
            G_con_12 = contravariant_metric[1, 2, i, j, element]
            G_con_22 = contravariant_metric[2, 2, i, j, element]
            for l in 1:2, m in 1:2
                christoffel_symbols[l, m, 1, i, j, element] = G_con_11 * Gamma_1[l, m] +
                                                              G_con_12 * Gamma_2[l, m]
                christoffel_symbols[l, m, 2, i, j, element] = G_con_12 * Gamma_1[l, m] +
                                                              G_con_22 * Gamma_2[l, m]
            end
        end
    end
end

# Check that the metric identities ∂ⱼ(J Gⁱʲ) = - JGʲᵏΓⱼₖⁱ are satisfied
function check_metric_compatibility(solver::DG,
                                    elements::P4estElementContainerCovariant{2},
                                    i, j, element)
    (; basis) = solver
    (; derivative_matrix) = basis
    (; contravariant_metric, christoffel_symbols) = elements

    J = volume_element(elements, i, j, element)
    JG_con_11 = J * contravariant_metric[1, 1, i, j, element]
    JG_con_12 = J * contravariant_metric[1, 2, i, j, element]
    JG_con_21 = J * contravariant_metric[2, 1, i, j, element]
    JG_con_22 = J * contravariant_metric[2, 2, i, j, element]

    # Differentiate with respect to ξ¹
    d1_JG_con_11 = zero(eltype(contravariant_metric))
    d1_JG_con_21 = zero(eltype(contravariant_metric))

    for ii in eachnode(basis)
        J_ii = volume_element(elements, ii, j, element)
        d1_JG_con_11 = d1_JG_con_11 +
                       derivative_matrix[i, ii] *
                       J_ii * contravariant_metric[1, 1, ii, j, element]
        d1_JG_con_21 = d1_JG_con_21 +
                       derivative_matrix[i, ii] *
                       J_ii * contravariant_metric[2, 1, ii, j, element]
    end

    # Differentiate with respect to ξ²
    d2_JG_con_12 = zero(eltype(contravariant_metric))
    d2_JG_con_22 = zero(eltype(contravariant_metric))
    for jj in eachnode(basis)
        J_jj = volume_element(elements, i, jj, element)
        d2_JG_con_12 = d2_JG_con_12 +
                       derivative_matrix[j, jj] *
                       J_jj * contravariant_metric[1, 2, i, jj, element]
        d2_JG_con_22 = d2_JG_con_22 +
                       derivative_matrix[j, jj] *
                       J_jj * contravariant_metric[2, 2, i, jj, element]
    end

    LHS1 = d1_JG_con_11 + d2_JG_con_12
    LHS2 = d1_JG_con_21 + d2_JG_con_22

    RHS1 = -(JG_con_11 * christoffel_symbols[1, 1, 1, i, j, element] +
             JG_con_12 * christoffel_symbols[1, 2, 1, i, j, element] +
             JG_con_21 * christoffel_symbols[2, 1, 1, i, j, element] +
             JG_con_22 * christoffel_symbols[2, 2, 1, i, j, element])

    RHS2 = -(JG_con_11 * christoffel_symbols[1, 1, 2, i, j, element] +
             JG_con_12 * christoffel_symbols[1, 2, 2, i, j, element] +
             JG_con_21 * christoffel_symbols[2, 1, 2, i, j, element] +
             JG_con_22 * christoffel_symbols[2, 2, 2, i, j, element])

    return LHS1 - RHS1, LHS2 - RHS2
end
end # @muladd
