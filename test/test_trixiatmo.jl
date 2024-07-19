# TODO - This is a repetition from Trixi.jl

using Test: @test
using Trixi: Trixi, examples_dir, trixi_include
import TrixiAtmo

"""
    @trixiatmo_testset "name of the testset" #= code to test #=

Similar to `@testset`, but wraps the code inside a temporary module to avoid
namespace pollution. It also `include`s this file again to provide the
definition of `@test_trixi_include`. Moreover, it records the execution time
of the testset similarly to [`timed_testset`](@ref).
"""
macro trixiatmo_testset(name, expr)
    @assert name isa String
    # TODO: `@eval` is evil
    # We would like to use
    #   mod = gensym(name)
    #   ...
    #   module $mod
    # to create new module names for every test set. However, this is not
    # compatible with the dirty hack using `@eval` to get the mapping when
    # loading structured, curvilinear meshes. Thus, we need to use a plain
    # module name here.
    quote
        local time_start = time_ns()
        @eval module TrixiAtmoTestModule
        using Test
        using TrixiAtmo
        include(@__FILE__)
        # We define `EXAMPLES_DIR` in (nearly) all test modules and use it to
        # get the path to the elixirs to be tested. However, that's not required
        # and we want to fail gracefully if it's not defined.
        try
            import ..EXAMPLES_DIR
            import ..TRIXI_EXAMPLES_DIR
        catch
            nothing
        end
        @testset $name $expr
        end
        local time_stop = time_ns()
        if TrixiAtmo.Trixi.mpi_isroot()
            flush(stdout)
            @info("Testset "*$name*" finished in "
                  *string(1.0e-9 * (time_stop - time_start))*" seconds.\n")
        end
        nothing
    end
end
