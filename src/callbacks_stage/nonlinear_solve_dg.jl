using Trixi: wrap_array, AbstractSemidiscretization#, TimerOutputs
#using TimerOutputs


@muladd begin

    
struct NonlinearSolver
    residual
    tolerance             ::Real
    variables_index_vector::Vector{Int}
end

#=
function NonlinearSolver(residual!, variables_index_vector; tolerance)
    NonlinearSolver(residual!, tolerance, variables_index_vector)
end
=#

function (limiter!::NonlinearSolver)(u_ode, integrator, semi::AbstractSemidiscretization, t)
    u = wrap_array(u_ode, semi)
    
    #@trixi_timeit timer() "nonlinear system solver" begin
        nonlinear_solve_dg_2d!(u, limiter!.residual, limiter!.tolerance, limiter!.variables_index_vector,
                               semi.equations, semi.solver, semi.cache)
    #end
end


include("nonlinear_solve_dg2d.jl")


end