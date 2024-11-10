# This method is called when a SemidiscretizationHyperbolic is constructed.
# It constructs the basic `cache` used throughout the simulation to compute
# the RHS etc.
function Trixi.create_cache(mesh::P4estMesh, equations::AbstractEquations, dg::DG, ::Any,
                            metric_terms, ::Type{uEltype}) where {uEltype <: Real}
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

    return cache
end

include("containers_2d_manifold_in_3d_cartesian.jl")
include("dg_2d_manifold_in_3d_cartesian.jl")

include("containers_2d_manifold_in_3d_covariant.jl")
include("dg_2d_manifold_in_3d_covariant.jl")
