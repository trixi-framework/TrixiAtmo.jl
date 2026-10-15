# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Abstract base type for the Rosenbrock-W methods with `N` stages, implemented in the
# transformed form of Hairer and Wanner (1996), Section IV.7: every stage solves one linear
# system with the matrix `I / (gamma * dt) - J`, where `J` is the approximate Jacobian of
# the right-hand side passed to `Trixi.solve` (`jacobian.jl`), and no product with the
# Jacobian is needed.
abstract type AbstractRosenbrockAlgorithm{N} <: AbstractTimeIntegrationAlgorithm end

nstages(::AbstractRosenbrockAlgorithm{N}) where {N} = N

# This struct is needed to fake https://github.com/SciML/OrdinaryDiffEq.jl/blob/0c2048a502101647ac35faabd80da8a5645beac7/src/integrators/type.jl#L77
# This implements the interface components described at
# https://diffeq.sciml.ai/v6.8/basics/integrator/#Handing-Integrators-1
# which are used in Trixi.jl.
mutable struct RosenbrockIntegrator{RealT <: Real, uType, Params, Sol, F, Alg,
                                    Options, NSTAGES, Jacobian} <:
               AbstractTimeIntegrator
    u::uType
    du::uType
    u_tmp::uType
    # Stage vectors; `stages[i]` holds the i-th stage of the transformed form
    stages::NTuple{NSTAGES, uType}
    t::RealT
    tdir::RealT # DIRection of time integration, i.e., if one marches forward or backward in time
    dt::RealT # current time step
    dtcache::RealT # manually set time step
    iter::Int # current number of time steps (iteration)
    p::Params # will be the semidiscretization from Trixi.jl
    sol::Sol # faked
    const f::F # right-hand side of the semidiscretization
    const alg::Alg # AbstractRosenbrockAlgorithm
    const jacobian::Jacobian # approximate Jacobian of the right-hand side, see `jacobian.jl`
    opts::Options
    finalstep::Bool # added for convenience
    const dtchangeable::Bool
    const force_stepfail::Bool
end

function init(ode::ODEProblem, alg::AbstractRosenbrockAlgorithm;
              dt, callback::Union{CallbackSet, Nothing} = nothing, jacobian = nothing,
              kwargs...)
    u = copy(ode.u0)
    du = similar(u)
    u_tmp = similar(u)
    # The stage vectors must be initialized to zero, not just allocated: the stage
    # combination sweeps run over all stages and multiply the stages not computed yet
    # by a zero coefficient, and `0 * NaN` would poison the result.
    stages = ntuple(_ -> zero(u), Val(nstages(alg)))
    t = first(ode.tspan)
    t, dt = promote(t, dt)
    tdir = sign(ode.tspan[end] - ode.tspan[1])
    iter = 0
    integrator = RosenbrockIntegrator(u, du, u_tmp, stages, t, tdir, dt, dt, iter,
                                      ode.p, (prob = ode,), ode.f, alg, jacobian,
                                      TimeIntegratorOptions(callback, ode.tspan;
                                                            kwargs...),
                                      false, true, false)

    # Standard callbacks
    initialize_callbacks!(callback, integrator)

    return integrator
end

function step!(integrator::RosenbrockIntegrator)
    @unpack prob = integrator.sol
    @unpack alg = integrator
    t_end = last(prob.tspan)
    callbacks = integrator.opts.callback

    @assert !integrator.finalstep
    if isnan(integrator.dt)
        error("time step size `dt` is NaN")
    end

    limit_dt!(integrator, t_end)

    # The matrix `I / (gamma dt) - J` of the stages is assembled once per time step. It
    # follows `dt` if `limit_dt!` shortens the last step.
    gamma_dt = alg.gamma * integrator.dt
    @trixi_timeit timer() "update Jacobian" update_jacobian!(integrator.jacobian,
                                                             integrator.u,
                                                             integrator.p, gamma_dt)

    for stage in 1:nstages(alg)
        rosenbrock_stage!(integrator, stage, gamma_dt)
    end

    integrator.iter += 1
    integrator.t += integrator.dt

    @trixi_timeit timer() "Step-Callbacks" handle_callbacks!(callbacks, integrator)

    check_max_iter!(integrator)

    return nothing
end

# One stage of the transformed form (Hairer and Wanner, 1996, eq. (7.25)):
#   u_s = u_n + sum_j a_sj stages_j
#   (I / (gamma dt) - J) stages_s = f(u_s, t + c_s dt) + sum_j (c_sj / dt) stages_j.
# After the last stage, u_{n+1} = u_n + sum_j m_j stages_j is written to `integrator.u`.
# The stage vector is assembled in place in `stages[stage]`: `res` is only an alias, so no
# separate residual vector is needed. This is safe because the tableau is strictly lower
# triangular, so the stage combination never uses `stages[stage]` (the GPU kernel skips the
# load of stages with a zero coefficient, and the CPU loop reads `stages[j][i]` for all `j`
# before writing `res[i]`).
function rosenbrock_stage!(integrator::RosenbrockIntegrator, stage, gamma_dt)
    (; u, du, u_tmp, stages, p, t, dt, alg) = integrator
    backend = trixi_backend(u)
    N = nstages(alg)

    a_row = alg.a[stage, :]
    c_row = alg.c[stage, :] * inv(dt)
    res = stages[stage]

    # Single fused sweep: u_tmp = u + sum_j a_sj stages_j, res = sum_j (c_sj / dt) stages_j
    @trixi_timeit_ext backend timer() "stage combination" rosenbrock_combine!(backend,
                                                                              u_tmp,
                                                                              res, u,
                                                                              stages,
                                                                              a_row,
                                                                              c_row,
                                                                              Val(N))

    integrator.f(du, u_tmp, p, t + alg.c_stage[stage] * dt)

    @trixi_timeit_ext backend timer() "stage residual" rosenbrock_residual!(backend,
                                                                            res, du)

    @trixi_timeit_ext backend timer() "solve Jacobian" solve_jacobian!(res,
                                                                       integrator.jacobian,
                                                                       p, gamma_dt)

    if stage == N
        # Single sweep: u = u + sum_j m_j stages_j (reads each stage once)
        m = alg.m
        @trixi_timeit_ext backend timer() "stage finalization" rosenbrock_finalize!(backend,
                                                                                    u,
                                                                                    stages,
                                                                                    m,
                                                                                    Val(N))
    end

    return nothing
end

# ---------------------------------------------------------------------------------------
# Full-field sweeps. These run over the whole solution vector and are the only parts of the
# integrator that touch `u` element by element, so they are the only parts that need a GPU
# counterpart. `trixi_backend(u)` returns `nothing` for a plain `Array` (threaded CPU loop)
# and a `KernelAbstractions.Backend` otherwise.
#
# The stage coefficients `a` and `c` are strictly lower triangular, so the sweeps always
# run over all `N` stages: the entries beyond the current stage are zero. This keeps the
# loop bound a compile-time constant, which lets it unroll on the GPU.
# ---------------------------------------------------------------------------------------

# u = u_n + sum_j a_j stages[j] and res = sum_j c_j stages[j], in one sweep
function rosenbrock_combine!(backend::Nothing, u, res, u_n, stages, a_row, c_row,
                             ::Val{N}) where {N}
    @threaded for i in eachindex(u)
        acc_u = u_n[i]
        acc_res = zero(eltype(res))
        for j in 1:N
            stage_ji = stages[j][i]
            acc_u = acc_u + a_row[j] * stage_ji
            acc_res = acc_res + c_row[j] * stage_ji
        end
        u[i] = acc_u
        res[i] = acc_res
    end
    return nothing
end

function rosenbrock_combine!(backend::Backend, u, res, u_n, stages, a_row, c_row,
                             ::Val{N}) where {N}
    kernel! = rosenbrock_combine_KAkernel!(backend)
    kernel!(u, res, u_n, stages, a_row, c_row, Val(N), ndrange = length(u))
    return nothing
end

# res = res + du
function rosenbrock_residual!(backend::Nothing, res, du)
    @threaded for i in eachindex(res)
        res[i] = res[i] + du[i]
    end
    return nothing
end

function rosenbrock_residual!(backend::Backend, res, du)
    res .+= du
    return nothing
end

# u = u + sum_j m_j stages[j], the final stage combination (in place)
function rosenbrock_finalize!(backend::Nothing, u, stages, m, ::Val{N}) where {N}
    @threaded for i in eachindex(u)
        acc = u[i]
        for j in 1:N
            acc = acc + m[j] * stages[j][i]
        end
        u[i] = acc
    end
    return nothing
end

function rosenbrock_finalize!(backend::Backend, u, stages, m, ::Val{N}) where {N}
    kernel! = rosenbrock_finalize_KAkernel!(backend)
    kernel!(u, stages, m, Val(N), ndrange = length(u))
    return nothing
end
end # @muladd

# sum_{j=1}^{J} coeffs[j] * stages[j][i], unrolled at compile time. The stages after the
# current one carry a zero coefficient; skipping the load rather than multiplying by zero
# saves the memory traffic, and the branch is uniform across the whole grid.
@inline function fold_stages(stages, coeffs, i, ::Val{J}) where {J}
    if J == 0
        return zero(eltype(coeffs))
    else
        previous = fold_stages(stages, coeffs, i, Val(J - 1))
        coeff = coeffs[J]
        return iszero(coeff) ? previous : coeff * (@inbounds stages[J][i]) + previous
    end
end

@kernel function rosenbrock_combine_KAkernel!(u, res, u_n, stages, a_row, c_row,
                                              ::Val{N}) where {N}
    i = @index(Global)
    @inbounds u[i] = u_n[i] + fold_stages(stages, a_row, i, Val(N))
    @inbounds res[i] = fold_stages(stages, c_row, i, Val(N))
end

@kernel function rosenbrock_finalize_KAkernel!(u, stages, m, ::Val{N}) where {N}
    i = @index(Global)
    @inbounds u[i] = u[i] + fold_stages(stages, m, i, Val(N))
end

include("tableaux.jl")
