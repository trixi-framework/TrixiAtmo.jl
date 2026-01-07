using Trixi: wrap_array, AbstractSemidiscretization, TimerOutputs, @trixi_timeit, timer

@muladd begin
#! format: noindent

"""
    NonlinearSolveDG

Newton method, which can be called in every stage via callbacks.

# Parameters
- `residual::Function`: function evaluating the residual
- `jacobian::Function`: function evaluating the Jacobian of the residual
- `variables_index_vector::Vector{Int64}`: vector of indices of entries of the solution vector the Newton method operates on
- `tolerance::Real`: tolerance for termination of the Newton method
- `max_iterations::Int64`: maximal number of iterations of the Newton method
"""
struct NonlinearSolveDG{RealT <: Real, ResidualFunctionT, JacobianFunctionT}
    residual               :: ResidualFunctionT
    jacobian               :: JacobianFunctionT
    variables_index_vector :: Vector{Int64}
    tolerance              :: RealT
    max_iterations         :: Int64

    function NonlinearSolveDG(residual, jacobian, variables_index_vector;
                              tolerance = 1e-9, max_iterations = 20)
        new{typeof(tolerance), typeof(residual), typeof(jacobian)}(residual, jacobian,
                                                                   variables_index_vector,
                                                                   tolerance,
                                                                   max_iterations)
    end
end

function (limiter!::NonlinearSolveDG)(u_ode, integrator,
                                      semi::AbstractSemidiscretization, t)
    u = wrap_array(u_ode, semi)

    @trixi_timeit timer() "nonlinear system solver" begin
        nonlinear_solve_dg_2d!(u, limiter!.residual, limiter!.jacobian,
                               limiter!.variables_index_vector,
                               limiter!.tolerance, limiter!.max_iterations,
                               semi.equations, semi.solver, semi.cache, semi.mesh)
    end
end
end

include("nonlinear_solve_dg2d.jl")
