@muladd begin
#! format: noindent

# New p4est element container that allows the use of a PtrArray for the 
# contravariant_vectors, and stores the auxiliary_variables as an additional field
mutable struct P4estElementContainerPtrArray{NDIMS, RealT <: Real, uEltype <: Real,
                                             NDIMSP1,
                                             NDIMSP2, NDIMSP3,
                                             ContravariantVectors <:
                                             AbstractArray{RealT,
                                                           NDIMSP3}} <:
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
    return size(elements.node_coordinates, ndims(elements) + 2)
end
@inline Base.ndims(::P4estElementContainerPtrArray{NDIMS}) where {NDIMS} = NDIMS
@inline function Base.eltype(::P4estElementContainerPtrArray{NDIMS,
                                                             RealT, uEltype}) where {NDIMS,
                                                                                     RealT,
                                                                                     uEltype}
    return uEltype
end

# Extract contravariant vector Ja^i (i = index) as SVector. This function dispatches on the 
# type of contravariant_vectors, specializing for NDIMS = 2 and NDIMS_AMBIENT = 3 by using 
# the fact that the second type parameter of PtrArray is NDIMS + 3, and the fourth type 
# parameter of PtrArray is Tuple{StaticInt{NDIMS_AMBIENT}, Vararg{IntT, NDIMS + 2}}.
@inline function Trixi.get_contravariant_vector(index,
                                                contravariant_vectors::PtrArray{RealT,
                                                                                5,
                                                                                <:Any,
                                                                                Tuple{Trixi.StaticInt{3},
                                                                                      IntT,
                                                                                      IntT,
                                                                                      IntT,
                                                                                      IntT}},
                                                indices...) where {RealT, IntT}
    return SVector(ntuple(@inline(dim->contravariant_vectors[dim, index, indices...]),
                          3))
end

# Create element container and initialize element data.
# This function dispatches on the dimensions of the mesh and the equation (AbstractEquations{3})
function Trixi.init_elements(mesh::Union{P4estMesh{2, 3, RealT},
                                         T8codeMesh{2}},
                             equations::AbstractEquations{3},
                             basis,
                             metric_terms,
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

    ContravariantVectors = typeof(contravariant_vectors)
    elements = P4estElementContainerPtrArray{NDIMS, RealT, uEltype, NDIMS + 1,
                                             NDIMS + 2,
                                             NDIMS + 3,
                                             ContravariantVectors}(node_coordinates,
                                                                   jacobian_matrix,
                                                                   contravariant_vectors,
                                                                   inverse_jacobian,
                                                                   surface_flux_values,
                                                                   _node_coordinates,
                                                                   _jacobian_matrix,
                                                                   _contravariant_vectors,
                                                                   _inverse_jacobian,
                                                                   _surface_flux_values)

    init_elements_2d_manifold_in_3d!(elements, mesh, basis, metric_terms)
    return elements
end

# This assumes a sphere, although the radius can be arbitrary, and general mesh topologies are supported.
function init_elements_2d_manifold_in_3d!(elements,
                                          mesh::P4estMesh{2, 3},
                                          basis::LobattoLegendreBasis,
                                          metric_terms)
    (; node_coordinates, jacobian_matrix,
     contravariant_vectors, inverse_jacobian) = elements

    # The standard calc_node_coordinates! can be used, since Trixi.jl now dispatches on
    # P4estMesh{NDIMS, NDIMS_AMBIENT}.
    Trixi.calc_node_coordinates!(node_coordinates, mesh, basis)

    for element in 1:Trixi.ncells(mesh)

        # Compute Jacobian matrix as usual
        Trixi.calc_jacobian_matrix!(jacobian_matrix, element, node_coordinates, basis)

        # Compute contravariant vectors
        calc_contravariant_vectors_2d_shell!(contravariant_vectors,
                                             element,
                                             jacobian_matrix,
                                             node_coordinates,
                                             basis,
                                             metric_terms)

        # Compute the inverse Jacobian as the norm of the cross product of the covariant vectors
        for j in eachnode(basis), i in eachnode(basis)
            inverse_jacobian[i, j, element] = 1 /
                                              norm(cross(jacobian_matrix[:, 1, i, j,
                                                                         element],
                                                         jacobian_matrix[:, 2, i, j,
                                                                         element]))
        end
    end

    return nothing
end

# Calculate contravariant vectors, multiplied by the Jacobian determinant J of the transformation mapping,
# using eq (12) of :
#   Giraldo, F. X. (2001). A spectral element shallow water model on spherical geodesic grids. 
#   International Journal for Numerical Methods in Fluids, 35(8), 869-901. 
#   https://doi.org/10.1002/1097-0363(20010430)35:8<869::AID-FLD116>3.0.CO;2-S
# This is equivalent to the cross-product form, but we end up with three contravariant vectors
# because there are three space dimensions.
# This only works for a sphere
function calc_contravariant_vectors_2d_shell!(contravariant_vectors::AbstractArray{<:Any,
                                                                                   5},
                                              element,
                                              jacobian_matrix, node_coordinates,
                                              basis::LobattoLegendreBasis,
                                              metric_terms::MetricTermsCrossProduct)
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

# Calculate contravariant vectors, multiplied by the Jacobian determinant J of the transformation mapping,
# using the invariant curl form for a cubed sphere.
# In the cubed sphere, the node coordinates are first mapped with polynomials in ξ and η, and then we add 
# a linear ζ dependency, i.e.:
#    xₗ(ξ, η, ζ) = ζ * xₗ(ξ, η) with ζ = 1 at the sphere surface.
# 
# As a result, the covariant vectors with respect to ζ are xₗ_ζ = xₗ
#
# We compute the metric terms in two steps
# 1. Compute the 3D contravariant vectors, Jaⁱ=J*grad(ξⁱ), using the curl invariant form and xₗ_ζ = xₗ
# 2. Convert the 3D mapping Jacobian determinant J:=a_k.(a_i x a_j) to 2D J_2D=||a_i x a_j||
## References
# - Kopriva, D. A. (2006). Metric identities and the discontinuous spectral element method on curvilinear meshes. Journal of #   Scientific Computing 26, 301-327. https://doi.org/10.1007/s10915-005-9070-8
# - Vinokur, M. and Yee, H. C. (2001). Extension of efficient low dissipation high order schemes for 3-D curvilinear moving # #   grids. In Caughey, D. A., and Hafez, M. M. (eds.), Frontiers of Computational Fluid Dynamics 2002, World Scientific, Singapore, pp. 129–164. https://doi.org/10.1142/9789812810793_0008

function calc_contravariant_vectors_2d_shell!(contravariant_vectors::AbstractArray{<:Any,
                                                                                   5},
                                              element,
                                              jacobian_matrix, node_coordinates,
                                              basis::LobattoLegendreBasis,
                                              metric_terms::MetricTermsInvariantCurl)
    @unpack derivative_matrix = basis

    # 1. Compute the 3D contravariant vectors, Jaⁱₙ=J*grad(ξ), using the curl invariant form and xₗ_ζ = xₗ
    # The general form is
    # Jaⁱₙ = 0.5 * ( ∇ × (Xₘ ∇ Xₗ - Xₗ ∇ Xₘ) )ᵢ  where (n, m, l) cyclic and ∇ = (∂/∂ξ, ∂/∂η, ∂/∂ζ)ᵀ

    for n in 1:3
        # (n, m, l) cyclic
        m = (n % 3) + 1
        l = ((n + 1) % 3) + 1

        # Calculate Ja¹ₙ = 0.5 * [ (Xₘ Xₗ_ζ - Xₗ Xₘ_ζ)_η - (Xₘ Xₗ_η - Xₗ Xₘ_η)_ζ ]
        # For each of these, the first and second summand are computed in separate loops
        # for performance reasons.

        # First summand 0.5 * (Xₘ Xₗ_ζ - Xₗ Xₘ_ζ)_η
        @turbo for j in eachnode(basis), i in eachnode(basis)
            result = zero(eltype(contravariant_vectors))

            for ii in eachnode(basis)
                # Multiply derivative_matrix to j-dimension to differentiate wrt η
                result += 0.5 * derivative_matrix[j, ii] *
                          (node_coordinates[m, i, ii, element] *
                           node_coordinates[l, i, ii, element] -
                           node_coordinates[l, i, ii, element] *
                           node_coordinates[m, i, ii, element])
            end

            contravariant_vectors[n, 1, i, j, element] = result
        end

        # Second summand -0.5 * (Xₘ Xₗ_η - Xₗ Xₘ_η)_ζ
        @turbo for j in eachnode(basis), i in eachnode(basis)
            # Due to the spherical-shell mapping, xₗ(ξ, η, ζ) = ζ * xₗ(ξ, η), we have
            # 0.5 d/dζ (ζ^2 F(ξ, η)) = ζ F(ξ, η)

            result = (node_coordinates[m, i, j, element] *
                      jacobian_matrix[l, 2, i, j, element] -
                      node_coordinates[l, i, j, element] *
                      jacobian_matrix[m, 2, i, j, element])

            contravariant_vectors[n, 1, i, j, element] -= result
        end

        # Calculate Ja²ₙ = 0.5 * [ (Xₘ Xₗ_ξ - Xₗ Xₘ_ξ)_ζ - (Xₘ Xₗ_ζ - Xₗ Xₘ_ζ)_ξ ]

        # First summand 0.5 * (Xₘ Xₗ_ξ - Xₗ Xₘ_ξ)_ζ
        @turbo for j in eachnode(basis), i in eachnode(basis)
            # Due to the spherical-shell mapping, xₗ(ξ, η, ζ) = ζ * xₗ(ξ, η), we have
            # 0.5 d/dζ (ζ^2 F(ξ, η)) = ζ F(ξ, η)

            result = (node_coordinates[m, i, j, element] *
                      jacobian_matrix[l, 1, i, j, element] -
                      node_coordinates[l, i, j, element] *
                      jacobian_matrix[m, 1, i, j, element])

            contravariant_vectors[n, 2, i, j, element] = result
        end

        # Second summand -0.5 * (Xₘ Xₗ_ζ - Xₗ Xₘ_ζ)_ξ
        @turbo for j in eachnode(basis), i in eachnode(basis)
            result = zero(eltype(contravariant_vectors))

            for ii in eachnode(basis)
                # Multiply derivative_matrix to i-dimension to differentiate wrt ξ
                result += 0.5 * derivative_matrix[i, ii] *
                          (node_coordinates[m, ii, j, element] *
                           node_coordinates[l, ii, j, element] -
                           node_coordinates[l, ii, j, element] *
                           node_coordinates[m, ii, j, element])
            end

            contravariant_vectors[n, 2, i, j, element] -= result
        end

        # Calculate Ja³ₙ = 0.5 * [ (Xₘ Xₗ_η - Xₗ Xₘ_η)_ξ - (Xₘ Xₗ_ξ - Xₗ Xₘ_ξ)_η ]

        # First summand 0.5 * (Xₘ Xₗ_η - Xₗ Xₘ_η)_ξ
        @turbo for j in eachnode(basis), i in eachnode(basis)
            result = zero(eltype(contravariant_vectors))

            for ii in eachnode(basis)
                # Multiply derivative_matrix to i-dimension to differentiate wrt ξ
                result += 0.5 * derivative_matrix[i, ii] *
                          (node_coordinates[m, ii, j, element] *
                           jacobian_matrix[l, 2, ii, j, element] -
                           node_coordinates[l, ii, j, element] *
                           jacobian_matrix[m, 2, ii, j, element])
            end

            contravariant_vectors[n, 3, i, j, element] = result
        end

        # Second summand -0.5 * (Xₘ Xₗ_ξ - Xₗ Xₘ_ξ)_η
        @turbo for j in eachnode(basis), i in eachnode(basis)
            result = zero(eltype(contravariant_vectors))

            for ii in eachnode(basis)
                # Multiply derivative_matrix to j-dimension to differentiate wrt η
                result += 0.5 * derivative_matrix[j, ii] *
                          (node_coordinates[m, i, ii, element] *
                           jacobian_matrix[l, 1, i, ii, element] -
                           node_coordinates[l, i, ii, element] *
                           jacobian_matrix[m, 1, i, ii, element])
            end

            contravariant_vectors[n, 3, i, j, element] -= result
        end
    end

    # 2. Convert the 3D mapping Jacobian determinant J:=a_k.(a_i x a_j) to 2D J_2D=||a_i x a_j||
    for j in eachnode(basis), i in eachnode(basis)
        factor = 1 / norm(node_coordinates[:, i, j, element]) # = 1 / norm(jacobian_matrix[:, 3, i, j, element])
        for n in 1:3, m in 1:3
            contravariant_vectors[m, n, i, j, element] *= factor
        end
    end

    return contravariant_vectors
end
end # @muladd
