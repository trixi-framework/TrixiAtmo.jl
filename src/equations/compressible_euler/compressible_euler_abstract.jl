# Abstract type for CompressibleEulerAtmo equations
#
# Generic functions go here
# Files including concrete types are included at the end

abstract type AbstractCompressibleEulerAtmo{NDIMS, NVARS, NPASSIVE} <:
              AbstractCompressibleEulerEquations{NDIMS, NVARS} end


@muladd begin
#! format: noindent

@inline nvariables_active(::AbstractCompressibleEulerAtmo{NDIMS, NVARS, NPASSIVE}) where
    {NDIMS, NVARS, NPASSIVE} = NVARS - NPASSIVE

@inline nvariables_passive(::AbstractCompressibleEulerAtmo{NDIMS, NVARS, NPASSIVE}) where
    {NDIMS, NVARS, NPASSIVE} = NPASSIVE

# TODO total or dry air?
#@inline function density(u, ::AbstractCompressibleEulerAtmo{NDIMS}) where {NDIMS}
#    return u[1]
#end

# Momentum densities come first
@inline function velocity(u, equations::AbstractCompressibleEulerAtmo{NDIMS}) where {NDIMS}
    rho = density_total(u, equations)
    return SVector(ntuple(@inline(v->u[v]/rho), Val(NDIMS)))
end

# thermodynamic variable next (NDIMS+1)

# dry air density next
@inline function density_dry(u, ::AbstractCompressibleEulerAtmo{NDIMS}) where {NDIMS}
    return u[NDIMS+2]
end

# all species, besides passive, are included when computing total density 
@inline function density_total(u,
    ::AbstractCompressibleEulerAtmo{NDIMS, NVARS, NPASSIVE}) where {NDIMS, NVARS, NPASSIVE}
    return sum(u[SVector{NVARS-NPASSIVE-NDIMS-1}(NDIMS+2:NVARS-NPASSIVE)])
end

@inline function prim2density_total(prim,
    ::AbstractCompressibleEulerAtmo{NDIMS, NVARS, NPASSIVE}) where {NDIMS, NVARS, NPASSIVE}
    # fraction of dry air density
    r_d = 1 - sum(prim[SVector{NVARS-NPASSIVE-NDIMS-2}(NDIMS+3:NVARS-NPASSIVE)])
    return prim[NDIMS+2] / r_d
end

# Calculate kinetic energy for a conservative state
# TODO: times rho ?!
#@inline function energy_kinetic(u, equations::AbstractCompressibleEulerAtmo)
#    rho = density_total(u, equations)
#    velocities = velocity(u, equations)
#    return 0.5f0 * rho * dot(velocities, velocities)
#end

# Transform initial conditions
# By convention all initial conditions are given in primitive variables
# in the form (rho, u, v, [w, ] p)
function transform_initial_condition(initial_condition, ::AbstractCompressibleEulerAtmo)
    function initial_condition_transformed(x, t,
        equations::AbstractCompressibleEulerAtmo{NDIMS, NVARS}) where
        {NDIMS, NVARS}
        prim = initial_condition(x, t, equations)
        prim_transformed = SVector(prim[SVector{NDIMS+1}(2:NDIMS+2)]...,
                                   prim[1],
                                   ntuple(i->0, Val(NVARS-NDIMS-2))...)
        return prim2cons(prim_transformed, equations)
    end
    return initial_condition_transformed
end

# Transform source terms
# By convention refers to conservative variables
# in the form (rho, rho u, rho v, [rho w, ] rho X) with X being the thermodynamic quantitiy
function transform_source_terms(source_terms, ::AbstractCompressibleEulerAtmo)
    function source_terms_transformed(u, x, t,
        equations::AbstractCompressibleEulerAtmo{NDIMS, NVARS}) where
        {NDIMS, NVARS}
        # rearrange to match with convention
        u_ref = SVector(u[NDIMS+2],
                        u[SVector{NDIMS+1}(1:NDIMS+1)]...) 
        source_ref = source_terms(u_ref, x, t, equations)
        # rearrange back
        return SVector(source_ref[SVector{NDIMS+1}(2:NDIMS+2)]...,
                       source_ref[1],
                       ntuple(i->0, Val(NVARS-NDIMS-2))...)
    end
    return source_terms_transformed
end

end # @muladd


# Include concrete types of equations
include("compressible_euler_atmo.jl")

# Include thermodynamic parts of equations
include("td_abstract.jl")
