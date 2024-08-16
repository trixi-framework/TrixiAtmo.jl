# The SemidiscretizationHyperbolic constructor has been modified to remove the assertion that checks 
# the compatibility between the mesh dimensionality and the equations' dimensionality. 
# Instead, we now directly dispatch using the specific mesh type (P4estMesh{2}) for 2D meshes and 
# AbstractEquations{3} for 3D equations. This change is necessary to support the Cartesian implementation 
# of a 2D manifold embedded in a 3D space.
function Trixi.SemidiscretizationHyperbolic(mesh::P4estMesh{2},
                                            equations::AbstractEquations{3},
                                            initial_condition,
                                            solver;
                                            source_terms = nothing,
                                            boundary_conditions = boundary_condition_periodic,
                                            # `RealT` is used as real type for node locations etc.
                                            # while `uEltype` is used as element type of solutions etc.
                                            RealT = real(solver), uEltype = RealT,
                                            initial_cache = NamedTuple())
    cache = (; Trixi.create_cache(mesh, equations, solver, RealT, uEltype)...,
             initial_cache...)
    _boundary_conditions = Trixi.digest_boundary_conditions(boundary_conditions, mesh,
                                                            solver,
                                                            cache)

    SemidiscretizationHyperbolic{typeof(mesh), typeof(equations),
                                 typeof(initial_condition),
                                 typeof(_boundary_conditions), typeof(source_terms),
                                 typeof(solver), typeof(cache)}(mesh, equations,
                                                                initial_condition,
                                                                _boundary_conditions,
                                                                source_terms, solver,
                                                                cache)
end
