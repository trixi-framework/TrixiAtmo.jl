# The SemidiscretizationHyperbolic constructor has been modified to remove the assertion
# that checks the compatibility between the mesh dimensionality and the equations' 
# dimensionality. Instead, we now directly dispatch using the specific mesh type 
# (P4estMesh{2}) for 2D meshes and AbstractEquations{3} for 3D equations.
# This change is necessary to support the implementation of a 2D manifold embedded in a 3D 
# space. We also pass in the keyword arguments "metric_terms" and "auxiliary_field", which 
# specify the approach used to compute the metric terms as well as any auxiliary fields 
# (e.g. bottom topography or a background state) that may be needed by the solver but are 
# not stored in the solution variables or elsewhere in the cache.
function Trixi.SemidiscretizationHyperbolic(mesh::P4estMesh{2},
                                            equations::AbstractEquations{3},
                                            initial_condition,
                                            solver;
                                            source_terms = nothing,
                                            boundary_conditions = boundary_condition_periodic,
                                            # `RealT` is used as real type for node locations etc.
                                            # while `uEltype` is used as element type of solutions etc.
                                            RealT = real(solver), uEltype = RealT,
                                            metric_terms = MetricTermsCrossProduct(),
                                            auxiliary_field = nothing)
    cache = (;
             Trixi.create_cache(mesh, equations, solver, RealT, metric_terms,
                                auxiliary_field, uEltype)...)
    _boundary_conditions = Trixi.digest_boundary_conditions(boundary_conditions, mesh,
                                                            solver,
                                                            cache)

    Trixi.check_periodicity_mesh_boundary_conditions(mesh, _boundary_conditions)

    performance_counter = Trixi.PerformanceCounter()

    SemidiscretizationHyperbolic{typeof(mesh), typeof(equations),
                                 typeof(initial_condition),
                                 typeof(_boundary_conditions), typeof(source_terms),
                                 typeof(solver), typeof(cache)}(mesh, equations,
                                                                initial_condition,
                                                                _boundary_conditions,
                                                                source_terms, solver,
                                                                cache,
                                                                performance_counter)
end

function Trixi.SemidiscretizationHyperbolic(mesh::DGMultiMesh,
                                            equations::AbstractCovariantEquations,
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

    performance_counter = Trixi.PerformanceCounter()

    SemidiscretizationHyperbolic{typeof(mesh), typeof(equations),
                                 typeof(initial_condition),
                                 typeof(_boundary_conditions), typeof(source_terms),
                                 typeof(solver), typeof(cache)}(mesh, equations,
                                                                initial_condition,
                                                                _boundary_conditions,
                                                                source_terms, solver,
                                                                cache,
                                                                performance_counter)
end

# Constructor for SemidiscretizationHyperbolic for the covariant form. Requires 
# compatibility between the mesh and equations (i.e. the same `NDIMS` and `NDIMS_AMBIENT`)
# and sets the default metric terms to MetricTermsCovariantSphere.
function Trixi.SemidiscretizationHyperbolic(mesh::P4estMesh{NDIMS, NDIMS_AMBIENT},
                                            equations::AbstractCovariantEquations{NDIMS,
                                                                                  NDIMS_AMBIENT},
                                            initial_condition,
                                            solver;
                                            source_terms = nothing,
                                            boundary_conditions = boundary_condition_periodic,
                                            # `RealT` is used as real type for node locations etc.
                                            # while `uEltype` is used as element type of solutions etc.
                                            RealT = real(solver), uEltype = RealT,
                                            metric_terms = MetricTermsCovariantSphere(),
                                            auxiliary_field = nothing) where {NDIMS,
                                                                              NDIMS_AMBIENT}
    cache = (;
             Trixi.create_cache(mesh, equations, solver, RealT, metric_terms,
                                auxiliary_field, uEltype)...)
    _boundary_conditions = Trixi.digest_boundary_conditions(boundary_conditions, mesh,
                                                            solver,
                                                            cache)

    Trixi.check_periodicity_mesh_boundary_conditions(mesh, _boundary_conditions)

    performance_counter = Trixi.PerformanceCounter()

    SemidiscretizationHyperbolic{typeof(mesh), typeof(equations),
                                 typeof(initial_condition),
                                 typeof(_boundary_conditions), typeof(source_terms),
                                 typeof(solver), typeof(cache)}(mesh, equations,
                                                                initial_condition,
                                                                _boundary_conditions,
                                                                source_terms, solver,
                                                                cache,
                                                                performance_counter)
end
