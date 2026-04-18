# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@inline function convert_to_solution_variables(u, solution_variables, cache,
                                               have_aux_node_vars::False,
                                               equations)
    # Reinterpret the solution array as an array of conservative variables,
    # compute the solution variables via broadcasting, and reinterpret the
    # result as a plain array of floating point numbers
    return Array(reinterpret(eltype(u),
                             solution_variables.(reinterpret(SVector{nvariables(equations),
                                                                     eltype(u)}, u),
                                                 Ref(equations))))
end

@inline function convert_to_solution_variables(u, solution_variables, cache,
                                               have_aux_node_vars::True,
                                               equations)
    @unpack aux_node_vars = cache.aux_vars
    # Reinterpret the solution array as an array of conservative variables,
    # compute the solution variables via broadcasting, and reinterpret the
    # result as a plain array of floating point numbers
    return Array(reinterpret(eltype(u),
                             solution_variables.(reinterpret(SVector{nvariables(equations),
                                                                     eltype(u)}, u),
                                                 reinterpret(SVector{n_aux_node_vars(equations),
                                                                     eltype(aux_node_vars)},
                                                             aux_node_vars),
                                                 Ref(equations))))
end

function Trixi.save_solution_file(u, time, dt, timestep,
                            mesh::P4estMesh,
                            equations, dg::DG, cache,
                            solution_callback,
                            element_variables = Dict{Symbol, Any}(),
                            node_variables = Dict{Symbol, Any}();
                            system = "")
    @unpack output_directory, solution_variables, reference_solution = solution_callback

    # Filename based on current time step
    if isempty(system)
        filename = joinpath(output_directory, @sprintf("solution_%09d.h5", timestep))
    else
        filename = joinpath(output_directory,
                            @sprintf("solution_%s_%09d.h5", system, timestep))
    end

    # Convert to different set of variables if requested
    converted_variables, n_vars = convert_variables(u, equations, solution_variables)

    # Subtract reference solution
    if !isnothing(reference_solution)
        data = converted_variables - reference_solution
    else
        data = converted_variables
    end

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
            attributes(var)["name"] = varnames(solution_variables, equations)[v]
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

end