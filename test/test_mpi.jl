@testitem "MPI: elixir_gemein_bubble" setup=[Setup] tags=[:mpi] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "euler", "dry_air", "buoyancy",
                                 "elixir_gemein_bubble.jl"),
                        l2=[
                            9.104437114458848e-7,
                            1.8210536975490044e-5,
                            0.0004707887343135412,
                            0.0063400898518523935,
                            0.0,
                            0.0
                        ],
                        linf=[
                            1.0258941581242631e-5,
                            0.00020520634691933992,
                            0.006392782691233334,
                            0.07637640493339859,
                            0.0,
                            0.0
                        ],
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 1000)
end

@testitem "MPI: elixir_potential_temperature_vortex_shedding with Sleve" setup=[Setup] tags=[:mpi] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "euler", "dry_air", "global_circulation",
                                 "elixir_potential_temperature_vortex_shedding.jl"),
                        l2=[
                            0.00010574237327573876,
                            0.03613038289661806,
                            0.03612899533359696,
                            0.04776452699672134,
                            0.03182161771653344
                        ],
                        linf=[
                            0.0017419769830600007,
                            0.15502918734879956,
                            0.15502638751987746,
                            0.2577665509810418,
                            0.34429461410371687
                        ],
                        rtol=1e-9,
                        tspan=(0.0, 0.0001 * SECONDS_PER_DAY),
                        trees_per_cube_face=(3, 2), adapt_vertical_grid=Sleve(0.7, 0.8))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 1000)
    # Check partitioning (a total of 108 elements split into 4 partitions)
    local_nelems = nelements(solver, semi.cache)

    # Perform the parallel reductions and assign the unwrapped results back!
    nelems_min = Trixi.MPI.Allreduce!(Ref(local_nelems), Base.min, Trixi.mpi_comm())[]
    nelems_max = Trixi.MPI.Allreduce!(Ref(local_nelems), Base.max, Trixi.mpi_comm())[]

    @assert nelems_min == 26
    @assert nelems_max == 28
end
