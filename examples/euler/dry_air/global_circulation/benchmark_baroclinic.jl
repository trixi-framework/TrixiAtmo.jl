using Trixi
using Trixi: trixi_include, wrap_array, mesh_equations_solver_cache, trixi_backend,
             have_nonconservative_terms, nelements, nvariables
using TrixiAtmo
using KernelAbstractions
using Printf
using Dates

const POLYDEGS = (3, 4, 5, 6, )
const VERSIONS = (0, 1, 5)
const TREES_PER_CUBE_FACES = ((30, 8), (15, 4), (20, 13), (16, 8))

const STORAGE_TYPE = CuArray
const REAL_TYPE = Float64
const TURBO = true 

const NSAMPLES = 100        
const WARMUP_CALLS = 5     

const OUTDIR = joinpath(pkgdir(Trixi), "run", "results")

const ALL_CASES = (
    (name = "Baroclinic_instability_with_log",
     elixir = joinpath(@__DIR__,
                       "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
     flux = TrixiAtmo.flux_etec_waruszewski_etal_combined,
     nonconservative_flux = (flux_etec, flux_nonconservative_waruszewski_etal)),
    (name = "Baroclinic_instability_no_log",
     elixir = joinpath(@__DIR__,
                       "elixir_potential_temperature_baroclinic_instability_turbo.jl"),
     flux = flux_kennedy_gruber_souza_etal,
     nonconservative_flux = (flux_kennedy_gruber, flux_nonconservative_souza_etal)),
)

const CASES = ALL_CASES[[1]]

function setup_case(case, polydeg, trees_per_cube_face; turbo = false)
    volume_flux = turbo ? FluxTurbo(case.nonconservative_flux) : case.flux

    elixir_module = Module()
    trixi_include(elixir_module, case.elixir;
                  polydeg = polydeg,
                  trees_per_cube_face = trees_per_cube_face,
                  volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
                  storage_type = STORAGE_TYPE, real_type = REAL_TYPE,
                  callbacks = nothing, sol = nothing)

    ode = elixir_module.ode
    semi = ode.p
    u_ode = ode.u0
    du_ode = similar(ode.u0)
    u = wrap_array(u_ode, semi)
    du = wrap_array(du_ode, semi)
    fill!(du, zero(eltype(du)))
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)
    backend = trixi_backend(u)

    return (; u, du, u_ode, du_ode, mesh, equations, solver, cache, backend,
            n_elements = nelements(solver, cache),
            n_variables = nvariables(equations))
end

function make_runner(setup, version::Integer, turbo::Bool)
    version_val = Val(version)
    turbo_trait = turbo ? Trixi.True() : Trixi.False()
    nonconservative = have_nonconservative_terms(setup.equations)

    backend = setup.backend
    du = setup.du
    u = setup.u
    mesh = setup.mesh
    equations = setup.equations
    volume_integral = setup.solver.volume_integral
    solver = setup.solver
    cache = setup.cache

    return function ()
        Trixi.calc_volume_integral!(backend, du, u, mesh, nonconservative,
                                    equations, volume_integral, solver, cache,
                                    version_val, turbo_trait)
        return nothing
    end
end

function time_runner(f, backend)
    total = 0.0
    best = Inf
    gc_was_enabled = GC.enable(false)
    try
        for _ in 1:NSAMPLES
            t = @elapsed begin
                f()
                KernelAbstractions.synchronize(backend)
            end
            total += t
            best = min(best, t)
        end
    finally
        GC.enable(gc_was_enabled)
    end
    return (; mean = total / NSAMPLES, best)
end

function benchmark_versions(setup, versions, turbo)
    runners = Dict{Int, Any}()
    failed = Int[]

    for v in versions
        runner = make_runner(setup, v, turbo)
        try
            runner()
            KernelAbstractions.synchronize(setup.backend)
            runners[v] = runner
        catch err
            @warn "version $v failed to launch" exception = err
            push!(failed, v)
        end
    end

    isempty(runners) && return Dict{Int, NamedTuple}(), failed

    for (_, r) in runners, _ in 1:WARMUP_CALLS
        r()
    end
    KernelAbstractions.synchronize(setup.backend)

    results = Dict(v => time_runner(r, setup.backend) for (v, r) in runners)

    return results, failed
end

function main()
    rows = NamedTuple[]

    @show STORAGE_TYPE
    @show REAL_TYPE
    @show VERSIONS
    @show TURBO

    for case in CASES,
        trees_per_cube_face in TREES_PER_CUBE_FACES,
        polydeg in POLYDEGS

        GC.gc()
        case_label = TURBO ? case.name * " + FluxTurbo" : case.name

        setup = try
            setup_case(case, polydeg, trees_per_cube_face; turbo = TURBO)
        catch err
            message = sprint(showerror, err)
            if occursin("out of memory", lowercase(message))
                @printf("\n%s  polydeg %d  trees %s : does not fit\n",
                        case_label, polydeg, string(trees_per_cube_face))
            else
                @printf("\n%s  polydeg %d  trees %s : FAILED, not a memory limit\n",
                        case_label, polydeg, string(trees_per_cube_face))
                println(message)
            end
            continue
        end

        n_dofs = setup.n_elements * (polydeg + 1)^3
        @printf("\n=== %s   polydeg %d   trees %s   %d elements   %.2fM DOF ===\n",
                case_label, polydeg, string(trees_per_cube_face),
                setup.n_elements, n_dofs / 1e6)

        results, failed = benchmark_versions(setup, VERSIONS, TURBO)

        for v in failed
            @printf("%-10d  launch failed\n", v)
        end
        if isempty(results)
            setup = nothing
            continue
        end

        ref_version = minimum(keys(results))
        ref_time = results[ref_version].best
        ref_version == first(VERSIONS) ||
            @warn "reference is v$ref_version, not v$(first(VERSIONS))"

        println(rpad("version", 9), rpad("best [ns/dof]", 15),
                rpad("mean [ns/dof]", 15), "vs v$ref_version")

        for v in sort(collect(keys(results)))
            r = results[v]

            @printf("%-9d%-15.2f%-15.2f%s\n",
                    v,
                    r.best * 1e9 / n_dofs,
                    r.mean * 1e9 / n_dofs,
                    @sprintf("%.2fx", ref_time / r.best))

            push!(rows,
                  (case = case_label,
                   polydeg = polydeg,
                   trees_per_cube_face = string(trees_per_cube_face),
                   elements = setup.n_elements,
                   n_dofs = n_dofs,
                   version = v,
                   reference_version = ref_version,
                   time_best = r.best,
                   time_mean = r.mean,
                   ns_per_dof_best = r.best * 1e9 / n_dofs,
                   ns_per_dof_mean = r.mean * 1e9 / n_dofs,
                   speedup_vs_ref = ref_time / r.best))
        end

        setup = nothing
    end

    write_csv(rows)
    return rows
end

function write_csv(rows)
    isempty(rows) && (println("\nno rows, nothing written"); return)

    mkpath(OUTDIR)
    stamp = Dates.format(now(), "yyyymmdd-HHMMSS")
    filename = string("benchmark",
                      "_v", join(VERSIONS, "-"),
                      "_p", join(POLYDEGS, "-"),
                      "_turbo-", TURBO,
                      "_", CASES[1].name,
                      "_", stamp, ".csv")
    path = joinpath(OUTDIR, filename)

    open(path, "w") do io
        println(io, join(string.(keys(rows[1])), ","))
        for r in rows
            fields = map(values(r)) do value
                value isa AbstractFloat ? @sprintf("%.10g", value) : string(value)
            end
            println(io, join(fields, ","))
        end
    end
    println("\nwrote ", path)
    return path
end

main()
