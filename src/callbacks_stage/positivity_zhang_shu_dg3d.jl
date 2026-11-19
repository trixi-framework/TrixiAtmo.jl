# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# The positivity-preserving limiter of Zhang and Shu of Trixi.jl, i.e., the
# `PositivityPreservingLimiterZhangShu`, on the CPU and with the kernels of
# KernelAbstractions.jl, e.g., on GPUs, for the equations with turbulence variables, e.g.,
# for the `density_turbulent_kinetic_energy` of the `DeardorffTKE`
function Trixi.limiter_zhang_shu!(u, threshold::Real, variable, mesh::P4estMesh{3},
                                  equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D,
                                  dg::DGSEM, cache)
    backend = trixi_backend(u)
    (; weights) = dg.basis
    (; inverse_jacobian) = cache.elements
    threshold = convert(eltype(u), threshold)
    if backend === nothing
        Trixi.@threaded for element in Trixi.eachelement(dg, cache)
            @inbounds limiter_zhang_shu_per_element!(u, threshold, variable, equations,
                                                     dg,
                                                     weights, inverse_jacobian, element)
        end
    else
        Trixi.nelements(dg, cache) == 0 && return nothing
        kernel! = limiter_zhang_shu_KAkernel!(backend)
        kernel!(u, threshold, variable, equations, dg, weights, inverse_jacobian,
                ndrange = Trixi.nelements(dg, cache))
    end
    return nothing
end

@kernel inbounds=true function limiter_zhang_shu_KAkernel!(u, threshold, variable,
                                                           equations, dg, weights,
                                                           inverse_jacobian)
    element = @index(Global)
    limiter_zhang_shu_per_element!(u, threshold, variable, equations, dg, weights,
                                   inverse_jacobian, element)
end

# Scale the deviations from the mean value of the element such that the minimum of
# `variable` at the nodes is at least `threshold`, as `limiter_zhang_shu!` of Trixi.jl
Base.@propagate_inbounds function limiter_zhang_shu_per_element!(u, threshold, variable,
                                                                 equations, dg, weights,
                                                                 inverse_jacobian,
                                                                 element)
    value_min = typemax(eltype(u))
    for k in Trixi.eachnode(dg), j in Trixi.eachnode(dg), i in Trixi.eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, k, element)
        value_min = min(value_min, variable(u_node, equations))
    end

    if value_min < threshold
        u_mean = zero(Trixi.get_node_vars(u, equations, dg, 1, 1, 1, element))
        total_volume = zero(eltype(u))
        for k in Trixi.eachnode(dg), j in Trixi.eachnode(dg), i in Trixi.eachnode(dg)
            node_volume = weights[i] * weights[j] * weights[k] *
                          abs(inv(inverse_jacobian[i, j, k, element]))
            total_volume += node_volume
            u_mean += Trixi.get_node_vars(u, equations, dg, i, j, k, element) *
                      node_volume
        end
        u_mean = u_mean / total_volume

        # The coefficient is clipped at zero if the mean value is below `threshold`, e.g.,
        # by round-off errors
        value_mean = variable(u_mean, equations)
        theta = max(zero(value_mean),
                    (value_mean - threshold) / (value_mean - value_min))

        for k in Trixi.eachnode(dg), j in Trixi.eachnode(dg), i in Trixi.eachnode(dg)
            limit_node!(u, theta, u_mean, variable, equations, dg, i, j, k, element)
        end
    end
    return nothing
end

Base.@propagate_inbounds function limit_node!(u, theta, u_mean, variable, equations, dg,
                                              i, j, k, element)
    u_node = Trixi.get_node_vars(u, equations, dg, i, j, k, element)
    Trixi.set_node_vars!(u, theta * u_node + (1 - theta) * u_mean, equations, dg,
                         i, j, k, element)
    return nothing
end

# The density of the TKE is a conservative variable, so only this variable is limited and
# the other variables are unchanged
Base.@propagate_inbounds function limit_node!(u, theta, u_mean,
                                              variable::typeof(density_turbulent_kinetic_energy),
                                              equations, dg, i, j, k, element)
    v = nvariables(equations)
    u[v, i, j, k, element] = theta * u[v, i, j, k, element] + (1 - theta) * u_mean[v]
    return nothing
end
end # @muladd
