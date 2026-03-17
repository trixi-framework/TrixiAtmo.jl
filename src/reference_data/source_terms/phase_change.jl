# By convention refers to conservative variables as defined in the equations

@muladd begin
#! format: noindent

function source_terms_phase_change_generator(::MicrophysicsRelaxation)

    @inline function source_terms(u, x, t,
                                  equations::CompressibleEulerAtmo{NDIMS}) where
                                  {NDIMS}
        u0 = zero(eltype(u))
        Q_ph = phase_change_vapour_liquid(u, equations, equations.microphysics)

        return SVector(ntuple(i -> u0, Val(NDIMS+2))..., Q_ph, -Q_ph)
    end 
    return source_terms
end
end # @muladd
