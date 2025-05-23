# The SemidiscretizationHyperbolic constructor has been modified to remove the assertion
# that checks the compatibility between the mesh dimensionality and the equations' 
# dimensionality. Instead, we now directly dispatch using the specific mesh type 
# (P4estMesh{2}) for 2D meshes and AbstractEquations{3} for 3D equations or 
# AbstractCovariantEquations{2, 3} for 2D covariant equations in three-dimensional space. 
# This change is necessary to support the implementation of a 2D manifold embedded in a 3D 
# space. We also pass in the keyword arguments "metric_terms" and "auxiliary_field", which 
# specify the approach used to compute the metric terms as well as any auxiliary fields 
# (e.g. bottom topography or a background state) that may be needed by the solver but are 
# not stored in the solution variables or elsewhere in the cache.
function Trixi.SemidiscretizationHyperbolic(mesh::P4estMesh{2},
                                            equations::Union{AbstractEquations{3},
                                                             AbstractCovariantEquations{2,
                                                                                        3}},
                                            initial_condition,
                                            solver;
                                            source_terms = nothing,
                                            boundary_conditions = boundary_condition_periodic,
                                            # `RealT` is used as real type for node locations etc.
                                            # while `uEltype` is used as element type of solutions etc.
                                            RealT = real(solver), uEltype = RealT,
                                            initial_cache = NamedTuple(),
                                            metric_terms = MetricTermsCrossProduct(),
                                            auxiliary_field = nothing)
    cache = (;
             Trixi.create_cache(mesh, equations, solver, RealT, metric_terms,
                                auxiliary_field, uEltype)..., initial_cache...)
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
