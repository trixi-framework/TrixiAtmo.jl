using Trixi: wrap_array, TimerOutputs, @trixi_timeit, timer



@muladd begin
 
struct RainLimiterDG
end


function (limiter!::RainLimiterDG)(u_ode, integrator, semi::AbstractSemidiscretization, t)
    u = wrap_array(u_ode, semi)
    
    @trixi_timeit timer() "rain limiter" begin
        rain_limiter_dg2d!(u, semi.equations, semi.solver, semi.cache, semi.mesh)
    end
end


include("rain_limiter_dg2d.jl")

end