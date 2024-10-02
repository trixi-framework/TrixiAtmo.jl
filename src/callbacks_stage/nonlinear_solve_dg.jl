using Trixi: wrap_array, AbstractSemidiscretization, TimerOutputs, @trixi_timeit, timer
using NonlinearSolve

@muladd begin
struct NonlinearSolveDG
    nonlin_fct::NonlinearFunction{false}
    variables_index_vector
    tolerance::Float64
end

function NonlinearSolveDG(residual, jacobian, variables_index_vector, tolerance)
    nonlin_fct = NonlinearFunction{false}(residual; jac = jacobian)
    return NonlinearSolveDG(nonlin_fct, variables_index_vector, tolerance)
end

function (limiter!::NonlinearSolveDG)(u_ode, integrator,
                                        semi::AbstractSemidiscretization, t)
    u = wrap_array(u_ode, semi)
    _, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    @unpack nonlin_fct, variables_index_vector, tolerance = limiter!

    @trixi_timeit timer() "nonlinear system solver" begin
        nonlinear_solve_dg_2d!(u, nonlin_fct, variables_index_vector, tolerance, equations, solver, cache)
    end
end

include("nonlinear_solve_dg2d.jl")
end
