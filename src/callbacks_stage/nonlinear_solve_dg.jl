using Trixi: wrap_array, AbstractSemidiscretization, TimerOutputs, @trixi_timeit, timer



@muladd begin
 
struct NonlinearSolveDG
    residual              ::Function
    jacobian              ::Function
    variables_index_vector::Vector{Int}
    tolerance             ::Real
end


function (limiter!::NonlinearSolveDG)(u_ode, integrator, semi::AbstractSemidiscretization, t)
    u = wrap_array(u_ode, semi)
    
    @trixi_timeit timer() "nonlinear system solver" begin
        nonlinear_solve_dg_2d!(u, limiter!.residual, limiter!.jacobian, limiter!.variables_index_vector,
                               limiter!.tolerance, semi.equations, semi.solver, semi.cache, semi.mesh)
    end
end


include("nonlinear_solve_dg2d.jl")

end