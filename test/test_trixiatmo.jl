# Testing uses TrixiTest.jl
# Collect everything needed here and include from test_*.jl

using Test
using TrixiTest
using TrixiAtmo: examples_dir

EXAMPLES_DIR = examples_dir()

# Check whether we run CI in the cloud with Window or Mac, see also
# https://docs.github.com/en/actions/learn-github-actions/environment-variables
CI_ON_WINDOWS = (get(ENV, "GITHUB_ACTIONS", false) == "true") && Sys.iswindows()
CI_ON_MAC = (get(ENV, "GITHUB_ACTIONS", false) == "true") && Sys.isapple()

macro test_trixi_include(expr, args...)
    local add_to_additional_ignore_content = [
    # This is needed because we overwrite `Trixi.weak_form_kernel!`, e.g., here:
    # https://github.com/trixi-framework/TrixiAtmo.jl/blob/20e069e818f23ed4033ef65fe087175f07d235fa/examples/elixir_shallowwater_cartesian_advection_cubed_sphere.jl#L50
        r"WARNING: Method definition .* in module .* at .* overwritten .*.\n"
    ]
    args = append_to_kwargs(args, :additional_ignore_content,
                            add_to_additional_ignore_content)
    quote
        @test_trixi_include_base($(esc(expr)), $(args...))
    end
end
