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

@testitem "MPI: elixir_potential_temperature_baroclinic_instability" setup=[Setup] tags=[:mpi] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "euler", "dry_air", "global_circulation",
                                 "elixir_potential_temperature_baroclinic_instability.jl"),
                        # Expected errors are exactly the same as in the serial test!
                        l2=[
                            0.0030320086112184466,
                            0.05144411394023726,
                            0.05238807358575555,
                            0.012609443441337428,
                            0.8067575925333563
                        ],
                        linf=[
                            0.021751907036954288,
                            0.8297225783677513,
                            1.2589734934841315,
                            0.9103791727770297,
                            4.833051748224818
                        ],
                        tspan=(0.0, 0.01 * SECONDS_PER_DAY), trees_per_cube_face=(2, 2),
                        dt=20)
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities). `Trixi.solve!` calls the right-hand side inside the
    # timer section "main loop", so the first call outside of it creates new timer
    # sections, which allocates.
    TrixiAtmo.Trixi.rhs_hyperbolic!(similar(sol.u[end]), sol.u[end], semi, sol.t[end])
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 1000)
    # Check partitioning: the vertical Jacobian requires that the columns are kept together
    # on the ranks
    @test nelements(solver, semi.cache) % TrixiAtmo.nlayers(mesh) == 0
end
