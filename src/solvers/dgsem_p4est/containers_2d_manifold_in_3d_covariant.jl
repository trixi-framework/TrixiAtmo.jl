@muladd begin
#! format: noindent

# Compute the auxiliary metric terms for the covariant form, assuming that the
# domain is a spherical shell. We do not make any assumptions on the mesh topology.
function compute_auxiliary_variables!(elements, mesh::P4estMesh{2, 3},
                                      equations::AbstractCovariantEquations{2, 3}, dg)

    # The tree node coordinates are assumed to be on the spherical shell centred around the 
    # origin. Therefore, we can compute the radius once and use it throughout.
    radius = sqrt(mesh.tree_node_coordinates[1, 1, 1, 1]^2 +
                  mesh.tree_node_coordinates[2, 1, 1, 1]^2 +
                  mesh.tree_node_coordinates[3, 1, 1, 1]^2)

    for element in 1:Trixi.ncells(mesh)

        # Check that the degree of the mesh matches that of the solver
        @assert length(mesh.nodes) == nnodes(dg)

        # Extract the corner vertex positions from the P4estMesh
        v1 = Trixi.get_node_coords(mesh.tree_node_coordinates, equations, dg,
                                   1, 1, element)
        v2 = Trixi.get_node_coords(mesh.tree_node_coordinates, equations, dg,
                                   nnodes(dg), 1, element)
        v3 = Trixi.get_node_coords(mesh.tree_node_coordinates, equations, dg,
                                   nnodes(dg), nnodes(dg), element)
        v4 = Trixi.get_node_coords(mesh.tree_node_coordinates, equations, dg,
                                   1, nnodes(dg), element)

        for j in eachnode(dg), i in eachnode(dg)

            # get reference node coordinates
            xi1, xi2 = dg.basis.nodes[i], dg.basis.nodes[j]

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
            basis_covariant = 0.25f0 * scaling_factor *
                              SMatrix{2, 3}(-sinlon, 0, coslon, 0, 0, 1) *
                              SMatrix{3, 3}(a11, a21, a31, a12, a22, a32, a13, a23,
                                            a33) *
                              SMatrix{3, 4}(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3],
                                            v3[1], v3[2], v3[3], v4[1], v4[2], v4[3]) *
                              SMatrix{4, 2}(-1 + xi2, 1 - xi2, 1 + xi2, -1 - xi2,
                                            -1 + xi1, -1 - xi1, 1 + xi1, 1 - xi1)

            # Set variables in the cache
            elements.auxiliary_variables[1, i, j, element] = basis_covariant[1, 1]
            elements.auxiliary_variables[2, i, j, element] = basis_covariant[2, 1]
            elements.auxiliary_variables[3, i, j, element] = basis_covariant[1, 2]
            elements.auxiliary_variables[4, i, j, element] = basis_covariant[2, 2]
        end
    end

    return nothing
end

# Get Cartesian node positions, dispatching on the dimension of the equation system
@inline function Trixi.get_node_coords(x,
                                       ::AbstractCovariantEquations{NDIMS,
                                                                    NDIMS_AMBIENT},
                                       ::DG,
                                       indices...) where {NDIMS, NDIMS_AMBIENT}
    return SVector(ntuple(@inline(idx->x[idx, indices...]), NDIMS_AMBIENT))
end

# Return the auxiliary variables at a given volume node index
@inline function get_node_aux_vars(auxiliary_variables, equations, solver::DG,
                                   indices...)
    SVector(ntuple(@inline(v->auxiliary_variables[v, indices...]),
                   Val(nauxvars(equations))))
end

# Return the auxiliary variables at a given surface node index
@inline function get_surface_node_aux_vars(auxiliary_variables, equations, solver::DG,
                                           indices...)
    aux_vars_ll = SVector(ntuple(@inline(v->auxiliary_variables[1, v, indices...]),
                                 Val(nauxvars(equations))))
    aux_vars_rr = SVector(ntuple(@inline(v->auxiliary_variables[2, v, indices...]),
                                 Val(nauxvars(equations))))
    return aux_vars_ll, aux_vars_rr
end

# Interface container storing surface values of the auxiliary variables
mutable struct P4estInterfaceContainerVariableCoefficient{NDIMS, uEltype <: Real,
                                                          NDIMSP2} <:
               Trixi.AbstractContainer
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

@inline Trixi.ninterfaces(interfaces::P4estInterfaceContainerVariableCoefficient) = size(interfaces.neighbor_ids,
                                                                                         2)
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
    auxiliary_variables = Trixi.unsafe_wrap(Array, pointer(_auxiliary_variables),
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

# Initialize node_indices of interface container (this is copied from Trixi.jl)
@inline function Trixi.init_interface_node_indices!(interfaces::P4estInterfaceContainerVariableCoefficient{2},
                                                    faces, orientation, interface_id)
    # Iterate over primary and secondary element
    for side in 1:2
        # Align interface in positive coordinate direction of primary element.
        # For orientation == 1, the secondary element needs to be indexed backwards
        # relative to the interface.
        if side == 1 || orientation == 0
            # Forward indexing
            i = :i_forward
        else
            # Backward indexing
            i = :i_backward
        end

        if faces[side] == 0
            # Index face in negative x-direction
            interfaces.node_indices[side, interface_id] = (:begin, i)
        elseif faces[side] == 1
            # Index face in positive x-direction
            interfaces.node_indices[side, interface_id] = (:end, i)
        elseif faces[side] == 2
            # Index face in negative y-direction
            interfaces.node_indices[side, interface_id] = (i, :begin)
        else # faces[side] == 3
            # Index face in positive y-direction
            interfaces.node_indices[side, interface_id] = (i, :end)
        end
    end

    return interfaces
end
end # muladd
