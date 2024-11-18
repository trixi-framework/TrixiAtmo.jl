@inline function get_node_aux_vars(auxiliary_variables, equations, solver::DG, indices...)
    SVector(ntuple(@inline(v->auxiliary_variables[v, indices...]),
                   Val(nauxvars(equations))))
end

@inline function get_surface_node_aux_vars(auxiliary_variables, equations, solver::DG,
                                           indices...)
    # There is a cut-off at `n == 10` inside of the method
    # `ntuple(f::F, n::Integer) where F` in Base at ntuple.jl:17
    # in Julia `v1.5`, leading to type instabilities if
    # more than ten variables are used. That's why we use
    # `Val(...)` below.
    aux_vars_ll = SVector(ntuple(@inline(v->auxiliary_variables[1, v, indices...]),
                                 Val(nauxvars(equations))))
    aux_vars_rr = SVector(ntuple(@inline(v->auxiliary_variables[2, v, indices...]),
                                 Val(nauxvars(equations))))
    return aux_vars_ll, aux_vars_rr
end

# Compute the node positions and metric terms for the covariant form, assuming that the
# domain is a spherical shell. We do not make any assumptions on the mesh topology.
function compute_auxiliary_variables!(elements, mesh::P4estMesh{2, 3},
                                      equations::AbstractCovariantEquations{2, 3},
                                      dg::DG)
    (; inverse_jacobian, auxiliary_variables) = elements
    (; basis) = dg

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
            auxiliary_variables[1, i, j, element] = 1 / inverse_jacobian[i, j, element]
            auxiliary_variables[2:5, i, j, element] = A
        end
    end

    return nothing
end

mutable struct P4estInterfaceContainerVariableCoefficient{NDIMS, uEltype <: Real,
                                                          NDIMSP2} <: Trixi.AbstractContainer
    u::Array{uEltype, NDIMSP2}       # [primary/secondary, variable, i, j, interface]
    auxiliary_variables::Array{uEltype, NDIMSP2} # [primary/secondary, variable, i, j, interface]
    neighbor_ids::Matrix{Int}                   # [primary/secondary, interface]
    node_indices::Matrix{NTuple{NDIMS, Symbol}} # [primary/secondary, interface]

    # internal `resize!`able storage
    _u::Vector{uEltype}
    _auxiliary_variables::Vector{uEltype}
    _neighbor_ids::Vector{Int}
    _node_indices::Vector{NTuple{NDIMS, Symbol}}
end

@inline function Trixi.ninterfaces(interfaces::P4estInterfaceContainerVariableCoefficient)
    size(interfaces.neighbor_ids, 2)
end
@inline Base.ndims(::P4estInterfaceContainerVariableCoefficient{NDIMS}) where {NDIMS} = NDIMS

# Create interface container and initialize interface data.
function Trixi.init_interfaces(mesh::Union{P4estMesh, T8codeMesh},
                               equations::AbstractCovariantEquations, basis, elements)
    NDIMS = ndims(elements)
    uEltype = eltype(elements)

    # Initialize container
    n_interfaces = Trixi.count_required_surfaces(mesh).interfaces

    _u = Vector{uEltype}(undef,
                         2 * nvariables(equations) * nnodes(basis)^(NDIMS - 1) *
                         n_interfaces)
    u = Trixi.unsafe_wrap(Array, pointer(_u),
                          (2, nvariables(equations),
                           ntuple(_ -> nnodes(basis), NDIMS - 1)...,
                           n_interfaces))

    _auxiliary_variables = Vector{uEltype}(undef,
                                           2 * nauxvars(equations) *
                                           nnodes(basis)^(NDIMS - 1) *
                                           n_interfaces)
    auxiliary_variables = Trixi.unsafe_wrap(Array, pointer(_u),
                                            (2, nauxvars(equations),
                                             ntuple(_ -> nnodes(basis), NDIMS - 1)...,
                                             n_interfaces))

    _neighbor_ids = Vector{Int}(undef, 2 * n_interfaces)
    neighbor_ids = Trixi.unsafe_wrap(Array, pointer(_neighbor_ids), (2, n_interfaces))

    _node_indices = Vector{NTuple{NDIMS, Symbol}}(undef, 2 * n_interfaces)
    node_indices = Trixi.unsafe_wrap(Array, pointer(_node_indices), (2, n_interfaces))

    interfaces = P4estInterfaceContainerVariableCoefficient{NDIMS, uEltype, NDIMS + 2}(u,
                                                                                       auxiliary_variables,
                                                                                       neighbor_ids,
                                                                                       node_indices,
                                                                                       _u,
                                                                                       _auxiliary_variables,
                                                                                       _neighbor_ids,
                                                                                       _node_indices)

    Trixi.init_interfaces!(interfaces, mesh)

    return interfaces
end

function prolong2interfaces_auxiliary!(cache, mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                       equations, surface_integral, dg::DG)
    (; elements, interfaces) = cache
    index_range = eachnode(dg)

    Trixi.@threaded for interface in eachinterface(dg, cache)
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
            for v in axes(auxiliary_variables, 1)
                interfaces.auxiliary_variables[1, v, i, interface] = elements.auxiliary_variables[v,
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
            for v in axes(auxiliary_variables, 1)
                interfaces.auxiliary_variables[2, v, i, interface] = elements.auxiliary_variables[v,
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
