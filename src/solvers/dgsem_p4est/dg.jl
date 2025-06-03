# This method is called when a SemidiscretizationHyperbolic is constructed.
# It constructs the basic `cache` used throughout the simulation to compute
# the RHS etc, and is modified here to use custom metric terms as well as provide the 
# option to # use auxiliary variables. 
function Trixi.create_cache(mesh::P4estMesh, equations::AbstractEquations, dg::DG, ::Any,
                            metric_terms, auxiliary_field,
                            ::Type{uEltype}) where {uEltype <: Real}
    # Make sure to balance the `p4est` before creating any containers
    # in case someone has tampered with the `p4est` after creating the mesh
    Trixi.balance!(mesh)

    elements = Trixi.init_elements(mesh, equations, dg.basis, metric_terms, uEltype)
    interfaces = Trixi.init_interfaces(mesh, equations, dg.basis, elements)
    boundaries = Trixi.init_boundaries(mesh, equations, dg.basis, elements)
    mortars = Trixi.init_mortars(mesh, equations, dg.basis, elements)

    cache = (; elements, interfaces, boundaries, mortars)

    # Add specialized parts of the cache required to compute the volume integral etc.
    cache = (; cache...,
             Trixi.create_cache(mesh, equations, dg.volume_integral, dg, uEltype)...)
    cache = (; cache..., Trixi.create_cache(mesh, equations, dg.mortar, uEltype)...)

    # Add specialized parts of the cache for auxiliary node variables
    cache = (; cache...,
             create_cache_auxiliary(mesh, equations,
                                    have_aux_node_vars(equations),
                                    dg, elements, interfaces, metric_terms,
                                    auxiliary_field)...)
    return cache
end

# If there are auxiliary variables, initialize them
function create_cache_auxiliary(mesh, equations, have_aux_node_vars::True, dg, elements,
                                interfaces, metric_terms, auxiliary_field)
    auxiliary_variables = init_auxiliary_node_variables(mesh, equations, dg, elements,
                                                        interfaces, metric_terms,
                                                        auxiliary_field)
    return (; auxiliary_variables)
end

# Do nothing if there are no auxiliary variables
function create_cache_auxiliary(mesh, equations, have_aux_node_vars::False, dg, elements,
                                interfaces, metric_terms, auxiliary_field)
    return NamedTuple()
end

include("containers_2d_manifold_in_3d_cartesian.jl")
include("containers_2d_manifold_in_3d_covariant.jl")

include("dg_2d_manifold_in_3d_cartesian.jl")
include("dg_2d_manifold_in_3d_covariant.jl")
