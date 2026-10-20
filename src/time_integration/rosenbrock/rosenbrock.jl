# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Abstract base type for time integration schemes of Rosenbrock-W methods with `N` stages.
# They are implemented in the transformed form of Hairer and Wanner (1996), Section IV.7,
# see `RosenbrockTableau`: every stage solves one linear system with the matrix
# `I / (gamma * dt) - J`, where `J` is the approximate Jacobian of the right-hand side
# passed to `Trixi.solve` (see `jacobian.jl`), and no matrix-vector product with the Jacobian is needed.
# A concrete method provides `tableau(alg, RealT)`.
abstract type AbstractRosenbrockAlgorithm{N} <: AbstractTimeIntegrationAlgorithm end

nstages(::AbstractRosenbrockAlgorithm{N}) where {N} = N

# This struct is needed to fake https://github.com/SciML/OrdinaryDiffEq.jl/blob/0c2048a502101647ac35faabd80da8a5645beac7/src/integrators/type.jl#L77
# This implements the interface components described at
# https://diffeq.sciml.ai/v6.8/basics/integrator/#Handing-Integrators-1
# which are used in Trixi.jl.
mutable struct RosenbrockIntegrator{RealT <: Real, uType, Params, Sol, F, Alg,
                                    Tableau, Options, NSTAGES, Jacobian} <:
               AbstractTimeIntegrator
    u::uType
    du::uType
    u_tmp::uType
    stages::NTuple{NSTAGES, uType} # stage vectors of the transformed form
    t::RealT
    tdir::RealT # DIRection of time integration, i.e., if one marches forward or backward in time
    dt::RealT # current time step
    dtcache::RealT # manually set time step
    iter::Int # current number of time steps (iteration)
    p::Params # will be the semidiscretization from Trixi.jl
    sol::Sol # faked
    const f::F # right-hand side of the semidiscretization
    const alg::Alg # AbstractRosenbrockAlgorithm
    const tableau::Tableau # RosenbrockTableau of `alg` in the real type of `u`
    const jacobian::Jacobian # arrays of the approximate Jacobian, see `jacobian.jl`
    const update_interval::Int # number of time steps between two updates of the Jacobian
    dt_jacobian::RealT # time step of the current factorization of the Jacobian
    opts::Options
    finalstep::Bool # added for convenience
    const dtchangeable::Bool
    const force_stepfail::Bool
end

function init(ode::ODEProblem, alg::AbstractRosenbrockAlgorithm;
              dt, callback::Union{CallbackSet, Nothing} = nothing, jacobian = nothing,
              update_interval = 1, kwargs...)
    u = copy(ode.u0)
    du = similar(u)
    u_tmp = similar(u)
    # The stage vectors are initialized to zero since the stage combinations run over all
    # stages, multiplying the stages not computed yet by a zero coefficient
    stages = ntuple(_ -> zero(u), Val(nstages(alg)))
    tableau_ = tableau(alg, eltype(u))
    t = first(ode.tspan)
    t, dt = promote(t, dt)
    tdir = sign(ode.tspan[end] - ode.tspan[1])
    iter = 0
    # The arrays of the approximate Jacobian are created for the semidiscretization and the
    # array type of the solution
    jacobian_cache = create_jacobian_cache(jacobian, ode.p, u)
    integrator = RosenbrockIntegrator(u, du, u_tmp, stages, t, tdir, dt, dt, iter,
                                      ode.p, (prob = ode,), ode.f, alg, tableau_,
                                      jacobian_cache, update_interval, zero(dt),
                                      TimeIntegratorOptions(callback, ode.tspan;
                                                            kwargs...),
                                      false, true, false)

    # Standard callbacks
    initialize_callbacks!(callback, integrator)

    return integrator
end

function step!(integrator::RosenbrockIntegrator)
    @unpack prob = integrator.sol
    @unpack alg, tableau = integrator
    t_end = last(prob.tspan)
    callbacks = integrator.opts.callback

    @assert !integrator.finalstep
    if isnan(integrator.dt)
        error("time step size `dt` is NaN")
    end

    modify_dt_for_tstops!(integrator)

    limit_dt!(integrator, t_end)

    # The stages divide by `dt`, so they are skipped for a vanishing time step, e.g., for
    # `tspan = (0.0, 0.0)`, which leaves the solution unchanged
    if !iszero(integrator.dt)
        # The matrix `I / (gamma dt) - J` of the stages is assembled every
        # `update_interval` time steps and whenever the time step changes (e.g., if a time
        # stop or `limit_dt!` shortens the step), since it depends on `dt`
        gamma_dt = tableau.gamma * integrator.dt
        if integrator.iter % integrator.update_interval == 0 ||
           integrator.dt != integrator.dt_jacobian
            @trixi_timeit timer() "update Jacobian" update_jacobian!(integrator.jacobian,
                                                                     integrator.u,
                                                                     integrator.p,
                                                                     gamma_dt)
            integrator.dt_jacobian = integrator.dt
        end

        for stage in 1:nstages(alg)
            rosenbrock_stage!(integrator, stage, gamma_dt)
        end
    end

    integrator.iter += 1
    integrator.t += integrator.dt

    @trixi_timeit timer() "Step-Callbacks" handle_callbacks!(callbacks, integrator)

    check_max_iter!(integrator)

    return nothing
end

# used for AMR
function Base.resize!(integrator::RosenbrockIntegrator, new_size)
    # The arrays of the approximate Jacobian depend on the mesh
    if !isnothing(integrator.jacobian)
        throw(ArgumentError("resizing the integrator (e.g., for AMR) is not supported " *
                            "with an approximate Jacobian"))
    end

    resize!(integrator.u, new_size)
    resize!(integrator.du, new_size)
    resize!(integrator.u_tmp, new_size)
    # The stages must be zero, since the stage combinations multiply the stages not
    # computed yet by a zero coefficient
    for stage in integrator.stages
        resize!(stage, new_size)
        fill!(stage, zero(eltype(stage)))
    end

    return nothing
end

# One stage of the transformed form (Hairer and Wanner, 1996, Eq. (7.25)),
#   u_i = u_n + sum_j a_ij stages_j,
#   (I / (gamma dt) - J) stages_i = f(u_i, t + c_i dt) + sum_j (c_ij / dt) stages_j,
# and, after the last stage, u_{n+1} = u_n + sum_j m_j stages_j.
function rosenbrock_stage!(integrator::RosenbrockIntegrator, stage, gamma_dt)
    (; u, du, u_tmp, stages, p, t, dt, alg, tableau) = integrator
    backend = trixi_backend(u)
    N = nstages(alg)

    a_row = tableau.a[stage, :]
    c_row = tableau.c[stage, :] * inv(dt)
    res = stages[stage]

    @trixi_timeit_ext backend timer() "internal stage" rosenbrock_internal_stage!(backend,
                                                                                  u_tmp,
                                                                                  res,
                                                                                  u,
                                                                                  stages,
                                                                                  a_row,
                                                                                  c_row,
                                                                                  Val(N))

    integrator.f(du, u_tmp, p, t + tableau.c_stage[stage] * dt)

    @trixi_timeit_ext backend timer() "stage residual" rosenbrock_stage_residual!(backend,
                                                                                  res,
                                                                                  du)

    @trixi_timeit_ext backend timer() "solve Jacobian" solve_jacobian!(res,
                                                                       integrator.jacobian,
                                                                       p, gamma_dt)

    if stage == N
        @trixi_timeit_ext backend timer() "final stage" rosenbrock_final_stage!(backend,
                                                                                u,
                                                                                stages,
                                                                                tableau.m,
                                                                                Val(N))
    end

    return nothing
end

# u = u_n + sum_j a_j stages[j] and res = sum_j c_j stages[j] in one sweep
function rosenbrock_internal_stage!(backend::Nothing, u, res, u_n, stages, a_row,
                                    c_row, ::Val{N}) where {N}
    @threaded for i in eachindex(u)
        u_i = u_n[i]
        res_i = zero(eltype(res))
        for j in 1:N
            stage_ji = stages[j][i]
            u_i = u_i + a_row[j] * stage_ji
            res_i = res_i + c_row[j] * stage_ji
        end
        u[i] = u_i
        res[i] = res_i
    end
    return nothing
end

function rosenbrock_internal_stage!(backend::Backend, u, res, u_n, stages, a_row,
                                    c_row, ::Val{N}) where {N}
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        Trixi.check_axes(res, axes(u))
        Trixi.check_axes(u_n, axes(u))
        for stage in stages
            Trixi.check_axes(stage, axes(u))
        end
    end
    kernel! = rosenbrock_internal_stage_KAkernel!(backend)
    kernel!(u, res, u_n, stages, a_row, c_row, Val(N), ndrange = length(u))
    return nothing
end

# res = res + du
function rosenbrock_stage_residual!(backend::Nothing, res, du)
    @threaded for i in eachindex(res)
        res[i] = res[i] + du[i]
    end
    return nothing
end

function rosenbrock_stage_residual!(backend::Backend, res, du)
    res .+= du
    return nothing
end

# u = u + sum_j m_j stages[j]
function rosenbrock_final_stage!(backend::Nothing, u, stages, m,
                                 ::Val{N}) where {N}
    @threaded for i in eachindex(u)
        u_i = u[i]
        for j in 1:N
            u_i = u_i + m[j] * stages[j][i]
        end
        u[i] = u_i
    end
    return nothing
end

function rosenbrock_final_stage!(backend::Backend, u, stages, m,
                                 ::Val{N}) where {N}
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        for stage in stages
            Trixi.check_axes(stage, axes(u))
        end
    end
    kernel! = rosenbrock_final_stage_KAkernel!(backend)
    kernel!(u, stages, m, Val(N), ndrange = length(u))
    return nothing
end

@kernel inbounds=true function rosenbrock_internal_stage_KAkernel!(u, res, u_n, stages,
                                                                   a_row, c_row,
                                                                   ::Val{N}) where {N}
    i = @index(Global)
    u_i = u_n[i]
    res_i = zero(eltype(res))
    for j in 1:N
        if !iszero(a_row[j]) || !iszero(c_row[j])
            stage_ji = stages[j][i]
            u_i = u_i + a_row[j] * stage_ji
            res_i = res_i + c_row[j] * stage_ji
        end
    end
    u[i] = u_i
    res[i] = res_i
end

@kernel inbounds=true function rosenbrock_final_stage_KAkernel!(u, stages, m,
                                                                ::Val{N}) where {N}
    i = @index(Global)
    u_i = u[i]
    for j in 1:N
        u_i = u_i + m[j] * stages[j][i]
    end
    u[i] = u_i
end
end # @muladd

include("tableaux.jl")
