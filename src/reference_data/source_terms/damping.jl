# By convention refers to conservative variables
# in the form (rho, rho u, rho v, [rho w, ] rho X) with X being the thermodynamic quantit

@muladd begin
#! format: noindent

function source_terms_rayleigh_damping_generator(;
    v1_0, v2_0,
    bounds,
    alpha = 0.5f0)

    @inline function source_terms(u, x, t,
                                  ::AbstractCompressibleEulerAtmo{NDIMS, NVARS, NPASSIVE}) where
                                  {NDIMS, NVARS, NPASSIVE}
        rho_total = u[1] + sum(u[SVector{NVARS-NPASSIVE-NDIMS-2}(NDIMS+3:NVARS-NPASSIVE)])
        v1 = u[2] / rho_total
        v2 = u[3] / rho_total
        u0 = zero(eltype(u))
        tau_s = u0
        for bound in bounds
            if abs(x[bound[1]]) > abs(bound[2])
                tau_s += alpha *
                    sin(0.5f0 * (x[2] - bound[2]) / (bound[3] - bound[2]))^2
            end
        end
        return SVector(u0,
                       -tau_s * rho_total * (v1 - v1_0),
                       -tau_s * rho_total * (v2 - v2_0),
                       u0)
    end
    return source_terms
end
end # @muladd
