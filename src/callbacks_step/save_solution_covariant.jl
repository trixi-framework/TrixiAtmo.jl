# Saves all solution variables in addition to the relative vorticity
function Trixi.save_solution_file(u, time, dt, timestep,
                                  mesh::Union{Trixi.SerialTreeMesh, StructuredMesh,
                                              StructuredMeshView,
                                              UnstructuredMesh2D, Trixi.SerialP4estMesh,
                                              Trixi.SerialT8codeMesh},
                                  equations::AbstractCovariantEquations, dg::DG,
                                  cache,
                                  solution_callback,
                                  element_variables = Dict{Symbol, Any}(),
                                  node_variables = Dict{Symbol, Any}();
                                  system = "")
    (; output_directory, solution_variables) = solution_callback
    (; aux_node_vars) = cache.auxiliary_variables

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
                                                                         eltype(u)},
                                                                 u, aux_node_vars),
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
