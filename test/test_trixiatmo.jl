# Testing uses TrixiTest.jl
# Collect everything needed here and include from test_*.jl

using TrixiTest
using TrixiAtmo

EXAMPLES_DIR = TrixiAtmo.examples_dir()

macro test_trixi_include(expr, args...)
    local add_to_additional_ignore_content = []
    args = append_to_kwargs(args, :additional_ignore_content,
                            add_to_additional_ignore_content)
    quote
        @test_trixi_include_base($(esc(expr)), $(args...))
    end
end
