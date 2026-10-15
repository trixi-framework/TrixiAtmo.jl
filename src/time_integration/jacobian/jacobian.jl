# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Interface of the approximate Jacobians of the linearly implicit time integration methods

"""
    update_jacobian!(jacobian, u_ode, semi, gamma_dt)

Assemble (and factorize) the approximate Jacobian `jacobian` of the right-hand side of the
semidiscretization `semi` at the solution `u_ode` for the Rosenbrock-W methods, e.g.,
[`SSPRW32`](@ref), such that [`solve_jacobian!`](@ref) can solve linear systems with the
matrix `I / gamma_dt - J`, where `gamma_dt = gamma * dt`. For `jacobian = nothing`,
nothing is done.
"""
update_jacobian!(::Nothing, u_ode, semi, gamma_dt) = nothing

"""
    solve_jacobian!(b_ode, jacobian, semi, gamma_dt)

Solve the linear system `(I / gamma_dt - J) x = b_ode` in place, where `J` is the
approximate Jacobian assembled by [`update_jacobian!`](@ref) and `gamma_dt = gamma * dt`.
For `jacobian = nothing`, i.e., `J = 0`, this scales `b_ode` by `gamma_dt`.
"""
function solve_jacobian!(b_ode, ::Nothing, semi, gamma_dt)
    b_ode .*= gamma_dt
    return nothing
end
end # @muladd
