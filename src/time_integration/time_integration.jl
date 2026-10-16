# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# The time integrators of TrixiAtmo.jl follow the ones of Trixi.jl: the algorithms are
# subtypes of `Trixi.AbstractTimeIntegrationAlgorithm` and the integrators subtypes of
# `Trixi.AbstractTimeIntegrator`, such that everything generic is reused from Trixi.jl,
# i.e., `Trixi.solve`, `Trixi.solve!`, the handling of the callbacks
# (`initialize_callbacks!`, `handle_callbacks!`, `finalize_callbacks`), of the final time
# step (`limit_dt!`), and of the maximal number of iterations (`check_max_iter!`), as well
# as the solution `TimeIntegratorSolution`.
# The parts of the integrator interface of SciMLBase.jl that Trixi.jl defines for each
# of its integrators separately are defined here once for all time integrators of
# TrixiAtmo.jl; a method only specializes them if it needs something different.

# This struct is needed to fake https://github.com/SciML/OrdinaryDiffEq.jl/blob/0c2048a502101647ac35faabd80da8a5645beac7/src/integrators/type.jl#L1
# It is the same as the options of the time integrators of Trixi.jl, e.g.,
# `Trixi.SimpleIntegratorSSPOptions`.
mutable struct TimeIntegratorOptions{Callback, TStops}
    callback::Callback # callbacks; used in Trixi.jl
    const adaptive::Bool # whether the algorithm is adaptive; ignored
    dtmax::Float64 # ignored
    const maxiters::Int # maximal number of time steps
    tstops::TStops # tstops from https://diffeq.sciml.ai/v6.8/basics/common_solver_opts/#Output-Control-1; ignored
end

function TimeIntegratorOptions(callback, tspan; maxiters = typemax(Int), kwargs...)
    tstops_internal = BinaryHeap{eltype(tspan)}(FasterForward())
    # We add last(tspan) to make sure that the time integration stops at the end time
    push!(tstops_internal, last(tspan))
    # We add 2 * last(tspan) because add_tstop!(integrator, t) is only called by DiffEqCallbacks.jl if tstops contains a time that is larger than t
    # (https://github.com/SciML/DiffEqCallbacks.jl/blob/025dfe99029bd0f30a2e027582744528eb92cd24/src/iterative_and_periodic.jl#L92)
    push!(tstops_internal, 2 * last(tspan))
    return TimeIntegratorOptions{typeof(callback), typeof(tstops_internal)}(callback,
                                                                            false,
                                                                            Inf,
                                                                            maxiters,
                                                                            tstops_internal)
end

"""
    add_tstop!(integrator::AbstractTimeIntegrator, t)
Add a time stop during the time integration process.
This function is called after the periodic SaveSolutionCallback to specify the next stop to save the solution.
"""
function add_tstop!(integrator::AbstractTimeIntegrator, t)
    integrator.tdir * (t - integrator.t) < zero(integrator.t) &&
        error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
    # We need to remove the first entry of tstops when a new entry is added.
    # Otherwise, the simulation gets stuck at the previous tstop and dt is adjusted to zero.
    if length(integrator.opts.tstops) > 1
        pop!(integrator.opts.tstops)
    end
    return push!(integrator.opts.tstops, integrator.tdir * t)
end

# get a cache where the RHS can be stored
get_tmp_cache(integrator::AbstractTimeIntegrator) = (integrator.u_tmp,)

# some algorithms from DiffEq like FSAL-ones need to be informed when a callback has modified u
derivative_discontinuity!(integrator::AbstractTimeIntegrator, ::Bool) = false

# stop the time integration
function terminate!(integrator::AbstractTimeIntegrator)
    integrator.finalstep = true

    return nothing
end
end # @muladd

include("jacobian/jacobian.jl")
include("rosenbrock/rosenbrock.jl")
