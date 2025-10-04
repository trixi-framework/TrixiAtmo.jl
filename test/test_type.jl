module TestUnit

using TrixiAtmo
using Trixi
using Test

include("test_trixiatmo.jl")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "Test Type Stability" begin
    @timed_testset "Compressible Euler Potential Temperature 1D" begin
        for RealT in (Float32, Float64)
            equations = @inferred CompressibleEulerPotentialTemperatureEquations1D(RealT)
        end
    end
end

end
