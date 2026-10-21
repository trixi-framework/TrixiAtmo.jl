# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Abstract base type for time integration schemes of linear IMEX Runge-Kutta methods with
# `N` stages. They split the right-hand side `f(u)` into `f(u) - J u`, treated explicitly,
# and `J u`, treated implicitly, where `J` is the approximate Jacobian passed to
# `Trixi.solve` (see `jacobian.jl`). They are implemented in the form of Giraldo, Kelly,
# and Constantinescu (2013), Section 3.2: every stage solves one linear system with the
# matrix `I / (gamma * dt) - J`, see `LinearIMEXTableau`, and no matrix-vector product with
# the Jacobian is needed. A concrete method provides `tableau(alg, RealT)`.
abstract type AbstractLinearIMEXAlgorithm{N} <: AbstractTimeIntegrationAlgorithm end

nstages(::AbstractLinearIMEXAlgorithm{N}) where {N} = N

# This struct is needed to fake https://github.com/SciML/OrdinaryDiffEq.jl/blob/0c2048a502101647ac35faabd80da8a5645beac7/src/integrators/type.jl#L77
# This implements the interface components described at
# https://diffeq.sciml.ai/v6.8/basics/integrator/#Handing-Integrators-1
# which are used in Trixi.jl.
mutable struct LinearIMEXIntegrator{RealT <: Real, uType, Params, Sol, F, Alg,
                                    Tableau, Options, NSTAGES, Jacobian} <:
               AbstractTimeIntegrator
    u::uType
    du::uType
    u_tmp::uType
    stages::NTuple{NSTAGES, uType} # stage values u_i; the first one is `u`
    rhs_stages::NTuple{NSTAGES, uType} # right-hand sides f(u_i) of the stages
    t::RealT
    tdir::RealT # DIRection of time integration, i.e., if one marches forward or backward in time
    dt::RealT # current time step
    dtcache::RealT # manually set time step
    iter::Int # current number of time steps (iteration)
    p::Params # will be the semidiscretization from Trixi.jl
    sol::Sol # faked
    const f::F # right-hand side of the semidiscretization
    const alg::Alg # AbstractLinearIMEXAlgorithm
    const tableau::Tableau # LinearIMEXTableau of `alg` in the real type of `u`
    const jacobian::Jacobian # arrays of the approximate Jacobian, see `jacobian.jl`
    const update_interval::Int # number of time steps between two updates of the Jacobian
    dt_jacobian::RealT # time step of the current factorization of the Jacobian
    opts::Options
    finalstep::Bool # added for convenience
    const dtchangeable::Bool
    const force_stepfail::Bool
end

function init(ode::ODEProblem, alg::AbstractLinearIMEXAlgorithm;
              dt, callback::Union{CallbackSet, Nothing} = nothing, jacobian = nothing,
              update_interval = 1, kwargs...)
    u = copy(ode.u0)
    du = similar(u)
    u_tmp = similar(u)
    # The first stage is explicit, i.e., u_1 = u_n, so its stage value is `u` itself. The
    # other stage vectors are initialized to zero since the stage combinations run over all
    # stages, multiplying the stages not computed yet by a zero coefficient.
    stages = (u, ntuple(_ -> zero(u), Val(nstages(alg) - 1))...)
    rhs_stages = ntuple(_ -> zero(u), Val(nstages(alg)))
    tableau_ = tableau(alg, eltype(u))
    t = first(ode.tspan)
    t, dt = promote(t, dt)
    tdir = sign(ode.tspan[end] - ode.tspan[1])
    iter = 0
    # The arrays of the approximate Jacobian are created for the semidiscretization and the
    # array type of the solution
    jacobian_cache = create_jacobian_cache(jacobian, ode.p, u)
    integrator = LinearIMEXIntegrator(u, du, u_tmp, stages, rhs_stages, t, tdir, dt, dt,
                                      iter, ode.p, (prob = ode,), ode.f, alg, tableau_,
                                      jacobian_cache, update_interval, zero(dt),
                                      TimeIntegratorOptions(callback, ode.tspan;
                                                            kwargs...),
                                      false, true, false)

    # Standard callbacks
    initialize_callbacks!(callback, integrator)

    return integrator
end

function step!(integrator::LinearIMEXIntegrator)
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
            linear_imex_stage!(integrator, stage, gamma_dt)
        end
    end

    integrator.iter += 1
    integrator.t += integrator.dt

    @trixi_timeit timer() "Step-Callbacks" handle_callbacks!(callbacks, integrator)

    check_max_iter!(integrator)

    return nothing
end

# used for AMR
function Base.resize!(integrator::LinearIMEXIntegrator, new_size)
    # The arrays of the approximate Jacobian depend on the mesh
    if !isnothing(integrator.jacobian)
        throw(ArgumentError("resizing the integrator (e.g., for AMR) is not supported " *
                            "with an approximate Jacobian"))
    end

    resize!(integrator.u, new_size)
    resize!(integrator.du, new_size)
    resize!(integrator.u_tmp, new_size)
    # The stages must be zero, since the stage combinations multiply the stages not
    # computed yet by a zero coefficient. The first stage value is `u`, which is already
    # resized.
    for stage in Base.tail(integrator.stages)
        resize!(stage, new_size)
        fill!(stage, zero(eltype(stage)))
    end
    for stage in integrator.rhs_stages
        resize!(stage, new_size)
        fill!(stage, zero(eltype(stage)))
    end

    return nothing
end

# One stage of the linear IMEX method (Giraldo, Kelly, and Constantinescu, 2013, Eq. (3.7))
# with the auxiliary variable z_i = u_i + sum_j d_ij u_j,
#   (I / (gamma dt) - J) z_i = (u_n + dt sum_j a_ij f(u_j) + sum_j d_ij u_j) / (gamma dt),
#   u_i = z_i - sum_j d_ij u_j,
# and, after the last stage, u_{n+1} = u_n + dt sum_j b_j f(u_j) (Eq. (3.6b)).
function linear_imex_stage!(integrator::LinearIMEXIntegrator, stage, gamma_dt)
    (; u, stages, rhs_stages, p, t, dt, alg, tableau) = integrator
    backend = trixi_backend(u)
    N = nstages(alg)

    # The first stage is explicit, i.e., u_1 = u_n
    if stage > 1
        a_row = tableau.a[stage, :] * dt
        d_row = tableau.d[stage, :]
        res = stages[stage]

        @trixi_timeit_ext backend timer() "internal stage" linear_imex_internal_stage!(backend,
                                                                                       res,
                                                                                       u,
                                                                                       stages,
                                                                                       rhs_stages,
                                                                                       a_row,
                                                                                       d_row,
                                                                                       inv(gamma_dt),
                                                                                       Val(N))

        @trixi_timeit_ext backend timer() "solve Jacobian" solve_jacobian!(res,
                                                                           integrator.jacobian,
                                                                           p, gamma_dt)

        @trixi_timeit_ext backend timer() "stage value" linear_imex_stage_value!(backend,
                                                                                 res,
                                                                                 stages,
                                                                                 d_row,
                                                                                 Val(N))
    end

    integrator.f(rhs_stages[stage], stages[stage], p, t + tableau.c[stage] * dt)

    if stage == N
        @trixi_timeit_ext backend timer() "final stage" linear_imex_final_stage!(backend,
                                                                                 u,
                                                                                 rhs_stages,
                                                                                 tableau.b *
                                                                                 dt,
                                                                                 Val(N))
    end

    return nothing
end

# res = (u_n + sum_j a_j rhs_stages[j] + sum_j d_j stages[j]) * inv_gamma_dt
function linear_imex_internal_stage!(backend::Nothing, res, u_n, stages, rhs_stages,
                                     a_row, d_row, inv_gamma_dt, ::Val{N}) where {N}
    # Explicit bounds check, which allows us to assume inbounds access below
    @boundscheck begin
        Trixi.check_axes(u_n, axes(res))
        for stage in stages
            Trixi.check_axes(stage, axes(res))
        end
        for stage in rhs_stages
            Trixi.check_axes(stage, axes(res))
        end
    end
    @threaded for i in eachindex(res)
        @inbounds begin
            res_i = u_n[i]
            for j in 1:N
                # Skip the stages not computed yet
                if !iszero(a_row[j])
                    res_i = res_i + a_row[j] * rhs_stages[j][i]
                end
                if !iszero(d_row[j])
                    res_i = res_i + d_row[j] * stages[j][i]
                end
            end
            res[i] = res_i * inv_gamma_dt
        end
    end
    return nothing
end

function linear_imex_internal_stage!(backend::Backend, res, u_n, stages, rhs_stages,
                                     a_row, d_row, inv_gamma_dt, ::Val{N}) where {N}
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        Trixi.check_axes(u_n, axes(res))
        for stage in stages
            Trixi.check_axes(stage, axes(res))
        end
        for stage in rhs_stages
            Trixi.check_axes(stage, axes(res))
        end
    end
    kernel! = linear_imex_internal_stage_KAkernel!(backend)
    kernel!(res, u_n, stages, rhs_stages, a_row, d_row, inv_gamma_dt, Val(N),
            ndrange = length(res))
    return nothing
end

@kernel inbounds=true function linear_imex_internal_stage_KAkernel!(res, u_n, stages,
                                                                    rhs_stages, a_row,
                                                                    d_row, inv_gamma_dt,
                                                                    ::Val{N}) where {N}
    i = @index(Global)
    res_i = u_n[i]
    for j in 1:N
        # Skip the stages not computed yet
        if !iszero(a_row[j])
            res_i = res_i + a_row[j] * rhs_stages[j][i]
        end
        if !iszero(d_row[j])
            res_i = res_i + d_row[j] * stages[j][i]
        end
    end
    res[i] = res_i * inv_gamma_dt
end

# res = res - sum_j d_j stages[j]
function linear_imex_stage_value!(backend::Nothing, res, stages, d_row,
                                  ::Val{N}) where {N}
    # Explicit bounds check, which allows us to assume inbounds access below
    @boundscheck begin
        for stage in stages
            Trixi.check_axes(stage, axes(res))
        end
    end
    @threaded for i in eachindex(res)
        @inbounds begin
            res_i = res[i]
            for j in 1:N
                # Skip the stages not computed yet
                if !iszero(d_row[j])
                    res_i = res_i - d_row[j] * stages[j][i]
                end
            end
            res[i] = res_i
        end
    end
    return nothing
end

function linear_imex_stage_value!(backend::Backend, res, stages, d_row,
                                  ::Val{N}) where {N}
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        for stage in stages
            Trixi.check_axes(stage, axes(res))
        end
    end
    kernel! = linear_imex_stage_value_KAkernel!(backend)
    kernel!(res, stages, d_row, Val(N), ndrange = length(res))
    return nothing
end

@kernel inbounds=true function linear_imex_stage_value_KAkernel!(res, stages, d_row,
                                                                 ::Val{N}) where {N}
    i = @index(Global)
    res_i = res[i]
    for j in 1:N
        # Skip the stages not computed yet
        if !iszero(d_row[j])
            res_i = res_i - d_row[j] * stages[j][i]
        end
    end
    res[i] = res_i
end

# u = u + sum_j b_j rhs_stages[j]
function linear_imex_final_stage!(backend::Nothing, u, rhs_stages, b,
                                  ::Val{N}) where {N}
    # Explicit bounds check, which allows us to assume inbounds access below
    @boundscheck begin
        for stage in rhs_stages
            Trixi.check_axes(stage, axes(u))
        end
    end
    @threaded for i in eachindex(u)
        @inbounds begin
            u_i = u[i]
            for j in 1:N
                u_i = u_i + b[j] * rhs_stages[j][i]
            end
            u[i] = u_i
        end
    end
    return nothing
end

function linear_imex_final_stage!(backend::Backend, u, rhs_stages, b,
                                  ::Val{N}) where {N}
    # Explicit bounds check, which allows us to assume inbounds access in the kernel
    @boundscheck begin
        for stage in rhs_stages
            Trixi.check_axes(stage, axes(u))
        end
    end
    kernel! = linear_imex_final_stage_KAkernel!(backend)
    kernel!(u, rhs_stages, b, Val(N), ndrange = length(u))
    return nothing
end

@kernel inbounds=true function linear_imex_final_stage_KAkernel!(u, rhs_stages, b,
                                                                 ::Val{N}) where {N}
    i = @index(Global)
    u_i = u[i]
    for j in 1:N
        u_i = u_i + b[j] * rhs_stages[j][i]
    end
    u[i] = u_i
end
end # @muladd

include("tableaux.jl")
