# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Abstract supertype of the time integrators of TrixiAtmo.jl for dispatch
abstract type AbstractTimeIntegrator end

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

has_tstop(integrator::AbstractTimeIntegrator) = !isempty(integrator.opts.tstops)
first_tstop(integrator::AbstractTimeIntegrator) = first(integrator.opts.tstops)

"""
    modify_dt_for_tstops!(integrator::AbstractTimeIntegrator)
Modify the time-step size to match the time stops specified in integrator.opts.tstops.
To avoid adding OrdinaryDiffEq to the dependencies, this routine is a copy of
https://github.com/SciML/OrdinaryDiffEq.jl/blob/d76335281c540ee5a6d1bd8bb634713e004f62ee/src/integrators/integrator_utils.jl#L38-L54
as in Trixi.jl.
"""
function modify_dt_for_tstops!(integrator::AbstractTimeIntegrator)
    if has_tstop(integrator)
        tdir_t = integrator.tdir * integrator.t
        tdir_tstop = first_tstop(integrator)
        if integrator.opts.adaptive
            integrator.dt = integrator.tdir *
                            min(abs(integrator.dt), abs(tdir_tstop - tdir_t)) # step! to the end
        elseif iszero(integrator.dtcache) && integrator.dtchangeable
            integrator.dt = integrator.tdir * abs(tdir_tstop - tdir_t)
        elseif integrator.dtchangeable && !integrator.force_stepfail
            # always try to step! with dtcache, but lower if a tstop
            # however, if force_stepfail then don't set to dtcache, and no tstop worry
            integrator.dt = integrator.tdir *
                            min(abs(integrator.dtcache), abs(tdir_tstop - tdir_t)) # step! to the end
        end
    end

    return nothing
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

# get a cache where the RHS can be stored
get_du(integrator::AbstractTimeIntegrator) = integrator.du

# Forward integrator.stats.naccept to integrator.iter (see GitHub PR#771 of Trixi.jl)
function Base.getproperty(integrator::AbstractTimeIntegrator, field::Symbol)
    if field === :stats
        return (naccept = getfield(integrator, :iter),)
    end
    # general fallback
    return getfield(integrator, field)
end

# used by adaptive timestepping algorithms in DiffEq
@inline function set_proposed_dt!(integrator::AbstractTimeIntegrator, dt)
    (integrator.dt = dt; integrator.dtcache = dt)

    return nothing
end

# Required e.g. for `glm_speed_callback`
@inline function get_proposed_dt(integrator::AbstractTimeIntegrator)
    return integrator.dt
end

@inline function limit_dt!(integrator::AbstractTimeIntegrator, t_end)
    # if the next iteration would push the simulation beyond the end time, set dt accordingly
    if integrator.t + integrator.dt > t_end ||
       isapprox(integrator.t + integrator.dt, t_end)
        integrator.dt = t_end - integrator.t
        terminate!(integrator)
    end

    return nothing
end

function initialize_callbacks!(callbacks::Union{CallbackSet, Nothing},
                               integrator::AbstractTimeIntegrator)
    # initialize callbacks
    if callbacks isa CallbackSet
        foreach(callbacks.continuous_callbacks) do cb
            throw(ArgumentError("Continuous callbacks are unsupported."))
        end
        foreach(callbacks.discrete_callbacks) do cb
            return cb.initialize(cb, integrator.u, integrator.t, integrator)
        end
    end

    return nothing
end

function handle_callbacks!(callbacks::Union{CallbackSet, Nothing},
                           integrator::AbstractTimeIntegrator)
    # handle callbacks
    if callbacks isa CallbackSet
        foreach(callbacks.discrete_callbacks) do cb
            if cb.condition(integrator.u, integrator.t, integrator)
                cb.affect!(integrator)
            end
            return nothing
        end
    end

    return nothing
end

@inline function check_max_iter!(integrator::AbstractTimeIntegrator)
    # respect maximum number of iterations
    if integrator.iter >= integrator.opts.maxiters && !integrator.finalstep
        @warn "Interrupted. Larger maxiters is needed."
        terminate!(integrator)
    end

    return nothing
end

# Called by `Trixi.solve` after `init`
function Trixi.solve!(integrator::AbstractTimeIntegrator)
    (; prob) = integrator.sol

    integrator.finalstep = false

    @trixi_timeit timer() "main loop" while !integrator.finalstep
        step!(integrator)
    end

    finalize_callbacks(integrator)

    return TimeIntegratorSolution((first(prob.tspan), integrator.t),
                                  (prob.u0, integrator.u), prob)
end

# Interface required by DiffEqCallbacks.jl
function DiffEqBase.get_tstops(integrator::AbstractTimeIntegrator)
    return integrator.opts.tstops
end
function DiffEqBase.get_tstops_array(integrator::AbstractTimeIntegrator)
    return get_tstops(integrator).valtree
end
function DiffEqBase.get_tstops_max(integrator::AbstractTimeIntegrator)
    return maximum(get_tstops_array(integrator))
end

function finalize_callbacks(integrator::AbstractTimeIntegrator)
    callbacks = integrator.opts.callback

    if callbacks isa CallbackSet
        foreach(callbacks.discrete_callbacks) do cb
            return cb.finalize(cb, integrator.u, integrator.t, integrator)
        end
        foreach(callbacks.continuous_callbacks) do cb
            return cb.finalize(cb, integrator.u, integrator.t, integrator)
        end
    end

    return nothing
end
end # @muladd

include("jacobian/jacobian.jl")
include("rosenbrock/rosenbrock.jl")
include("linear_imex/linear_imex.jl")
