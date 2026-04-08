@muladd begin

mutable struct SaveSolutionCallbackRef{IntervalType, SolutionVariablesType,
                                    ReferenceSolutionType}
    interval_or_dt::IntervalType
    save_initial_solution::Bool
    save_final_solution::Bool
    output_directory::String
    solution_variables::SolutionVariablesType
    reference_solution::ReferenceSolutionType
end

# Convenience function to convert variables
function Trixi.save_solution_file(u, time, dt, timestep,
                            mesh:: P4estMesh{2},
                            equations::PerturbationEulerEquations2DAuxVars, dg::DG, cache,
                            solution_callback,
                            element_variables = Dict{Symbol, Any}(),
                            node_variables = Dict{Symbol, Any}();
                            system = "")
    #@unpack output_directory, solution_variables, reference_solution = solution_callback
    @unpack output_directory, solution_variables = solution_callback
    (; aux_node_vars) = cache.aux_vars
    # Filename based on current time step
    if isempty(system)
        filename = joinpath(output_directory, @sprintf("solution_%09d.h5", timestep))
    else
        filename = joinpath(output_directory,
                            @sprintf("solution_%s_%09d.h5", system, timestep))
    end

    # Convert to different set of variables if requested
    converted_variables, n_vars = convert_variables(u, aux_node_vars, equations, solution_variables)

    # Subtract reference solution
    #if !isnothing(reference_solution)
    #    data = converted_variables - reference_solution
    #else
    #    data = converted_variables
    #end
    data = converted_variables
    # Open file (clobber existing content)
    h5open(filename, "w") do file
        # Add context information as attributes
        attributes(file)["ndims"] = ndims(mesh)
        attributes(file)["equations"] = Trixi.get_name(equations)
        attributes(file)["polydeg"] = Trixi.polydeg(dg)
        attributes(file)["n_vars"] = n_vars
        attributes(file)["n_elements"] = Trixi.nelements(dg, cache)
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

function convert_variables(u, aux, equations::PerturbationEulerEquations2DAuxVars, solution_variables)
    if solution_variables === Trixi.cons2cons
        data = u
        n_vars = Trixi.nvariables(equations)
    else
        # Reinterpret the solution array as an array of conservative variables,
        # compute the solution variables via broadcasting, and reinterpret the
        # result as a plain array of floating point numbers
        data = Array(reinterpret(eltype(u),
                                 solution_variables.(reinterpret(SVector{Trixi.nvariables(equations),
                                                                         eltype(u)}, u), 
                                                     reinterpret(SVector{n_aux_node_vars(equations), eltype(u)}, aux),
                                                     Ref(equations))))

        # Find out variable count by looking at output from `solution_variables` function
        n_vars = size(data, 1)
    end

    return data, n_vars
end

# Convenience function to convert variables
function convert_variables(u, equations, solution_variables)
    if solution_variables === cons2cons
        data = u
        n_vars = nvariables(equations)
    else
        # Reinterpret the solution array as an array of conservative variables,
        # compute the solution variables via broadcasting, and reinterpret the
        # result as a plain array of floating point numbers
        data = Array(reinterpret(eltype(u),
                                 solution_variables.(reinterpret(SVector{nvariables(equations),
                                                                         eltype(u)}, u),
                                                     Ref(equations))))

        # Find out variable count by looking at output from `solution_variables` function
        n_vars = size(data, 1)
    end

    return data, n_vars
end

#for substracting a reference (the background state) from the solution
function SaveSolutionCallbackRef(; interval::Integer = 0,
                              dt = nothing,
                              save_initial_solution = true,
                              save_final_solution = true,
                              output_directory = "out",
                              solution_variables = cons2prim,
                              reference_solution = nothing)
    if !isnothing(dt) && interval > 0
        throw(ArgumentError("You can either set the number of steps between output (using `interval`) or the time between outputs (using `dt`) but not both simultaneously"))
    end

    # Expected most frequent behavior comes first
    if isnothing(dt)
        interval_or_dt = interval
    else # !isnothing(dt)
        interval_or_dt = dt
    end

    solution_callback = SaveSolutionCallbackRef(interval_or_dt,
                                             save_initial_solution, save_final_solution,
                                             output_directory, solution_variables,
                                             reference_solution)

    # Expected most frequent behavior comes first
    if isnothing(dt)
        # Save every `interval` (accepted) time steps
        # The first one is the condition, the second the affect!
        return DiscreteCallback(solution_callback, solution_callback,
                                save_positions = (false, false),
                                initialize = initialize_save_cb!)
    else
        # Add a `tstop` every `dt`, and save the final solution.
        return PeriodicCallback(solution_callback, dt,
                                save_positions = (false, false),
                                initialize = initialize_save_cb!,
                                final_affect = save_final_solution)
    end
end

# Compute a reference state to be subtracted from the solution at each write to file
function compute_reference_state(func, semi, solution_variables)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)

    reference_state = Trixi.compute_coefficients(func, 0.0, semi)
    reference_state_wrapped = copy(Trixi.wrap_array_native(reference_state,
                                                           mesh, equations, solver,
                                                           cache))

    data, _ = convert_variables(reference_state_wrapped, equations, solution_variables)
    return data
end


end 