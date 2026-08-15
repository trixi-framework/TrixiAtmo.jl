# By convention refers to conservative variables
# in the form (rho, rho u, rho v, [rho w, ] rho X) with X being the thermodynamic quantit

@muladd begin
#! format: noindent

function source_terms_rayleigh_damping_generator(;
                                                 v1_0, v2_0,
                                                 bounds,
                                                 alpha = 0.5f0)
    @inline function source_terms(u, x, t,
                                  equations::AbstractCompressibleEulerAtmo{2})
        rho_total = density_total(u, equations)
        v = vars_moment(u, equations) / rho_total
        tau_s = zero(eltype(u))
        for bound in bounds
            if abs(x[bound[1]]) > abs(bound[2])
                tau_s += alpha *
                         sinpi(0.5f0 * (x[bound[1]] - bound[2]) / (bound[3] - bound[2]))^2
            end
        end
        return SVector(-tau_s * rho_total * (v[1] - v1_0),
                       -tau_s * rho_total * (v[2] - v2_0))
    end
    return source_terms
end
end # @muladd
