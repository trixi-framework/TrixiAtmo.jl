using Trixi: wrap_array, AbstractSemidiscretization, TimerOutputs, @trixi_timeit, timer
using NonlinearSolve

@muladd begin
struct NonlinearSolveDG
    nonlin_fct
    variables_index_vector
    tolerance::Float64
end

function NonlinearSolveDG(residual, jacobian, variables_index_vector, tolerance)
    nonlin_fct = NonlinearFunction(residual; jac = jacobian)
    return NonlinearSolveDG(nonlin_fct, variables_index_vector, tolerance)
end

function (limiter!::NonlinearSolveDG)(u_ode, integrator,
                                        semi::AbstractSemidiscretization, t)
    u = wrap_array(u_ode, semi)

    @trixi_timeit timer() "nonlinear system solver" begin
        nonlinear_solve_dg_2d!(u, limiter!.nonlin_fct, limiter!.variables_index_vector,
                               limiter!.tolerance, semi.equations, semi.solver, semi.cache)
    end
end

include("nonlinear_solve_dg2d.jl")
end
