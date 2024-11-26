# Convert to another set of variables using a solution_variables function that depends on 
# auxiliary variables
function convert_variables(u, solution_variables, mesh::P4estMesh{2},
                           equations, dg, cache)
    (; aux_node_vars) = cache.auxiliary_variables

    # Extract the number of solution variables to be output 
    # (may be different than the number of conservative variables) 
    n_vars = length(solution_variables(Trixi.get_node_vars(u, equations, dg, 1, 1, 1),
                                       get_node_aux_vars(aux_node_vars, equations, dg, 1, 1,
                                                         1), equations))
    # Allocate storage for output
    data = Array{eltype(u)}(undef, n_vars, nnodes(dg), nnodes(dg), nelements(dg, cache))

    # Loop over all nodes and convert variables, passing in auxiliary variables
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
            a_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            u_node_converted = solution_variables(u_node, a_node, equations)
            for v in 1:n_vars
                data[v, i, j, element] = u_node_converted[v]
            end
        end
    end
    return data
end

# Version of save_solution_file that supports a solution_variables function that depends on 
# auxiliary variables
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
    else
        data = convert_variables(u, solution_variables, mesh, equations, dg, cache)
    end
    n_vars = size(data, 1)

    # Open file (clobber existing content)
    h5open(filename, "w") do file
        # Add context information as attributes
        attributes(file)["ndims"] = ndims(mesh)
        attributes(file)["equations"] = Trixi.get_name(equations)
        attributes(file)["polydeg"] = Trixi.polydeg(dg)
        attributes(file)["n_vars"] = n_vars
        attributes(file)["n_elements"] = nelements(dg, cache)
        attributes(file)["mesh_type"] = Trixi.get_name(mesh)
        attributes(file)["mesh_file"] = splitdir(mesh.current_filename)[2]
        attributes(file)["time"] = convert(Float64, time) # Ensure that `time` is written as a double precision scalar
        attributes(file)["dt"] = convert(Float64, dt) # Ensure that `dt` is written as a double precision scalar
        attributes(file)["timestep"] = timestep

        # Store each variable of the solution data
        for v in 1:n_vars
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
