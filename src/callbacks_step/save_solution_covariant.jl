@muladd begin
#! format: noindent

# Calculate the relative vorticity at all nodes in the mesh
function calc_relative_vorticity!(relative_vorticity, u, mesh::P4estMesh,
                                  equations::CovariantShallowWaterEquations2D, dg::DG,
                                  cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            relative_vorticity[i, j, element] = calc_relative_vorticity(u, equations,
                                                                        dg, cache, i, j,
                                                                        element)
        end
    end
end

# Calculate the relative vorticity at a given node using the collocation derivative
@inline function calc_relative_vorticity(u, equations::CovariantShallowWaterEquations2D,
                                         dg::DGSEM, cache, i, j, element)
    (; derivative_matrix) = dg.basis
    (; covariant_metric, inverse_jacobian) = cache.elements

    dv2dxi1 = dv1dxi2 = zero(eltype(u))

    for ii in eachnode(dg)
        h, hv_con_1, hv_con_2 = Trixi.get_node_vars(u, equations, dg, ii, j, element)
        hv_cov_2 = covariant_metric[2, 1, ii, j, element] * hv_con_1 +
                   covariant_metric[2, 2, ii, j, element] * hv_con_2
        dv2dxi1 = dv2dxi1 + derivative_matrix[i, ii] * hv_cov_2 / h
    end

    for jj in eachnode(dg)
        h, hv_con_1, hv_con_2 = Trixi.get_node_vars(u, equations, dg, i, jj, element)
        hv_cov_1 = covariant_metric[1, 1, i, jj, element] * hv_con_1 +
                   covariant_metric[1, 2, i, jj, element] * hv_con_2
        dv1dxi2 = dv1dxi2 + derivative_matrix[j, jj] * hv_cov_1 / h
    end

    # compute the relative vorticity
    return (dv2dxi1 - dv1dxi2) * inverse_jacobian[i, j, element]
end

# Saves all solution variables in addition to the relative vorticity
function Trixi.save_solution_file(u, time, dt, timestep,
                                  mesh::Union{Trixi.SerialTreeMesh, StructuredMesh,
                                              StructuredMeshView,
                                              UnstructuredMesh2D, Trixi.SerialP4estMesh,
                                              Trixi.SerialT8codeMesh},
                                  equations::CovariantShallowWaterEquations2D, dg::DG,
                                  cache,
                                  solution_callback,
                                  element_variables = Dict{Symbol, Any}(),
                                  node_variables = Dict{Symbol, Any}();
                                  system = "")
    (; output_directory, solution_variables) = solution_callback

    # Filename based on current time step
    if isempty(system)
        filename = joinpath(output_directory, @sprintf("solution_%09d.h5", timestep))
    else
        filename = joinpath(output_directory,
                            @sprintf("solution_%s_%09d.h5", system, timestep))
    end

    # Convert to different set of variables if requested
    if solution_variables === cons2cons
        data = u
        n_vars_cons = nvariables(equations)
    else
        # Reinterpret the solution array as an array of conservative variables,
        # compute the solution variables via broadcasting, and reinterpret the
        # result as a plain array of floating point numbers
        data = Array(reinterpret(eltype(u),
                                 solution_variables.(reinterpret(SVector{nvariables(equations),
                                                                         eltype(u)}, u),
                                                     Ref(equations))))

        # Find out variable count by looking at output from `solution_variables` function
        n_vars_cons = size(data, 1)
    end

    # Open file (clobber existing content)
    h5open(filename, "w") do file
        # Add context information as attributes
        attributes(file)["ndims"] = ndims(mesh)
        attributes(file)["equations"] = Trixi.get_name(equations)
        attributes(file)["polydeg"] = Trixi.polydeg(dg)
        attributes(file)["n_vars"] = n_vars_cons + 1
        attributes(file)["n_elements"] = nelements(dg, cache)
        attributes(file)["mesh_type"] = Trixi.get_name(mesh)
        attributes(file)["mesh_file"] = splitdir(mesh.current_filename)[2]
        attributes(file)["time"] = convert(Float64, time) # Ensure that `time` is written as a double precision scalar
        attributes(file)["dt"] = convert(Float64, dt) # Ensure that `dt` is written as a double precision scalar
        attributes(file)["timestep"] = timestep

        # Store each variable of the solution data
        for v in 1:n_vars_cons
            # Convert to 1D array
            file["variables_$v"] = vec(data[v, .., :])

            # Add variable name as attribute
            var = file["variables_$v"]
            attributes(var)["name"] = Trixi.varnames(solution_variables, equations)[v]
        end

        # Calculate relative vorticity
        relative_vorticity = similar(data[1, .., :])
        calc_relative_vorticity!(relative_vorticity, u, mesh, equations, dg, cache)
        nvarsp1 = n_vars_cons + 1
        file["variables_$nvarsp1"] = vec(relative_vorticity)
        attributes(file["variables_$nvarsp1"])["name"] = "relative_vorticity"

        # Store element variables
        for (v, (key, element_variable)) in enumerate(element_variables)
            # Add to file
            file["element_variables_$v"] = element_variable

            # Add variable name as attribute
            var = file["element_variables_$v"]
            attributes(var)["name"] = string(key)
        end

        # Store node variables
        for (v, (key, node_variable)) in enumerate(node_variables)
            # Add to file
            file["node_variables_$v"] = node_variable

            # Add variable name as attribute
            var = file["node_variables_$v"]
            attributes(var)["name"] = string(key)
        end
    end

    return filename
end
end # @muladd
