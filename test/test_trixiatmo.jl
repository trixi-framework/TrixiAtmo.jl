# Testing uses TrixiTest.jl
# Collect everything needed here and include from test_*.jl

using TrixiTest
using TrixiAtmo

EXAMPLES_DIR = TrixiAtmo.examples_dir()

macro test_trixi_include(expr, additional_ignore_content = [])
    quote
        add_to_additional_ignore_content = []
        append!($additional_ignore_content, add_to_additional_ignore_content)
        @test_trixi_include_base $(esc(expr)) $additional_ignore_content
    end
end
