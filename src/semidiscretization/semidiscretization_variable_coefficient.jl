"""
        function Trixi.SemidiscretizationHyperbolic(mesh::P4estMesh{2},
                                            equations::AbstractVariableCoefficientEquations{2,1},
                                            initial_condition,
                                            solver;
                                            source_terms = nothing, 
                                            boundary_conditions = boundary_condition_periodic,
                                            RealT = real(solver), uEltype = RealT,
                                            initial_cache = NamedTuple(),
                                            auxiliary_field = nothing)

`aux_field` is an optional function taking a coordinate vector `x` and the current
`equations` as arguments. It is used to fill an additional field `aux_vars` in the `cache`,
which will be available, e.g., in flux computations. The current `equations` need to set
`have_aux_node_vars` to `True()` and `n_aux_node_vars` to the number of auxiliary variables
per node.
"""


function Trixi.SemidiscretizationHyperbolic(mesh::P4estMesh{2},
                                            equations::AbstractVariableCoefficientEquations{},
                                            initial_condition,
                                            solver;
                                            source_terms = nothing, 
                                            boundary_conditions = boundary_condition_periodic,
                                            # `RealT` is used as real type for node locations etc.
                                            # while `uEltype` is used as element type of solutions etc.
                                            RealT = real(solver), uEltype = RealT,
                                            initial_cache = NamedTuple(),
                                            aux_field = nothing)
        cache = (; Trixi.create_cache(mesh, equations, solver, RealT, nothing, aux_field, uEltype)...,
                initial_cache...)
        _boundary_conditions = Trixi.digest_boundary_conditions(boundary_conditions, mesh, solver, cache)

        SemidiscretizationHyperbolic{typeof(mesh), typeof(equations),
                                     typeof(initial_condition),
                                     typeof(_boundary_conditions), typeof(source_terms),
                                     typeof(solver), typeof(cache)}(mesh, equations,
                                     initial_condition,
                                     _boundary_conditions,
                                     source_terms, solver,
                                     cache)
end