@muladd begin
#! format: noindent

# Convert to another set of variables using a solution_variables function that depends on 
# auxiliary variables
function convert_variables(u, solution_variables, mesh::P4estMesh{2},
                           equations::AbstractEquations{3}, dg, cache)
    (; contravariant_vectors) = cache.elements
    # Extract the number of solution variables to be output 
    # (may be different than the number of conservative variables) 
    n_vars = length(Trixi.varnames(solution_variables, equations))

    # Allocate storage for output
    data = Array{eltype(u)}(undef, n_vars, nnodes(dg), nnodes(dg), nelements(dg, cache))

    # Loop over all nodes and convert variables, passing in auxiliary variables
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            
            normal_vector_node = Trixi.get_contravariant_vector(3, contravariant_vectors, i, j, element)

            data_node = solution_variables(u, normal_vector_node, equations, dg, cache, i, j, element)

            for v in 1:n_vars
                data[v, i, j, element] = data_node[v]
            end
        end
    end
    return data
end

# If no auxiliary variables are passed into the conversion to spherical coordinates, do not 
# do any conversion.
@inline cons2prim_plus_vorticity(u, equations) = u

# # Specialized save_solution_file method that supports a solution_variables function which 
# # depends on auxiliary variables. The conversion must be defined as solution_variables(u, 
# # aux_vars, equations), and an additional method must be defined as solution_variables(u,
# # equations) = u, such that no conversion is done when auxiliary variables are not provided.
# function Trixi.save_solution_file(u_ode, t, dt, iter,
#                                   semi::SemidiscretizationHyperbolic{<:Trixi.AbstractMesh{2},
#                                                                      <:AbstractEquations{3}},
#                                   solution_callback,
#                                   element_variables = Dict{Symbol, Any}(),
#                                   node_variables = Dict{Symbol, Any}();
#                                   system = "")
#     mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
#     u = Trixi.wrap_array_native(u_ode, semi)

#     # Perform the variable conversion at each node
#     data = convert_variables(u, solution_callback.solution_variables, mesh, equations,
#                              solver, cache)

#     solution_callback2 = copy(solution_callback)
#     solution_callback2
#     # Call the existing Trixi.save_solution_file, which will use solution_variables(u, 
#     # equations). Since the variables are already converted above, we therefore require 
#     # solution_variables(u, equations) = u.
#     Trixi.save_solution_file(data, t, dt, iter, mesh, equations, solver, cache,
#                              solution_callback, element_variables,
#                              node_variables, system = system)
# end


# Calculate the relative vorticity at a given node using the collocation derivative
@inline function cons2prim_plus_vorticity(u, normal_vector, equations::ShallowWaterEquations3D,
                                          dg::DGSEM, cache, i, j, element)
    (; derivative_matrix) = dg.basis
    (; contravariant_vectors, inverse_jacobian) = cache.elements

    # Compute gradients in reference space
    dv1dxi1 = dv1dxi2 = zero(eltype(u))
    dv2dxi1 = dv2dxi2 = zero(eltype(u))
    dv3dxi1 = dv3dxi2 = zero(eltype(u))
    for ii in eachnode(dg)
        h, hv_1, hv_2, hv_3, _ = Trixi.get_node_vars(u, equations, dg, ii, j, element)
        dv1dxi1 += derivative_matrix[i, ii] * hv_1 / h
        dv2dxi1 += derivative_matrix[i, ii] * hv_2 / h
        dv3dxi1 += derivative_matrix[i, ii] * hv_3 / h
    end

    for jj in eachnode(dg)
        h, hv_1, hv_2, hv_3, _ = Trixi.get_node_vars(u, equations, dg, i, jj, element)
        dv1dxi2 += derivative_matrix[j, jj] * hv_1 / h
        dv2dxi2 += derivative_matrix[j, jj] * hv_2 / h
        dv3dxi2 += derivative_matrix[j, jj] * hv_3 / h
    end

    # Transform gradients to Cartesian space
    Ja11, Ja12, Ja13 = Trixi.get_contravariant_vector(1, contravariant_vectors, i, j, element)
    Ja21, Ja22, Ja23 = Trixi.get_contravariant_vector(2, contravariant_vectors, i, j, element)
    
    dv1dy = (Ja12 * dv1dxi1 + Ja22 * dv1dxi2) * inverse_jacobian[i, j, element]
    dv1dz = (Ja13 * dv1dxi1 + Ja23 * dv1dxi2) * inverse_jacobian[i, j, element]
    dv2dx = (Ja11 * dv2dxi1 + Ja21 * dv2dxi2) * inverse_jacobian[i, j, element]
    dv2dz = (Ja13 * dv2dxi1 + Ja23 * dv2dxi2) * inverse_jacobian[i, j, element]
    dv3dx = (Ja11 * dv3dxi1 + Ja21 * dv3dxi2) * inverse_jacobian[i, j, element]
    dv3dy = (Ja12 * dv3dxi1 + Ja22 * dv3dxi2) * inverse_jacobian[i, j, element]

    # compute the vorticity and project onto normal vector
    vorticity = ((dv3dy - dv2dz) * normal_vector[1] + 
                 (dv1dz - dv3dx) * normal_vector[2] + 
                 (dv2dx - dv1dy) * normal_vector[3]) / norm(normal_vector)

    u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)

    return SVector(cons2prim(u_node, equations)..., vorticity)
end

Trixi.varnames(::typeof(cons2prim_plus_vorticity), equations::ShallowWaterEquations3D) = (varnames(cons2prim, equations)..., "vorticity")

function Trixi.save_solution_file(u, time, dt, timestep,
                            mesh::P4estMesh{2},
                            equations::AbstractEquations{3}, dg::DG, cache,
                            solution_callback,
                            element_variables = Dict{Symbol, Any}(),
                            node_variables = Dict{Symbol, Any}();
                            system = "")
    @unpack output_directory, solution_variables = solution_callback

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
        n_vars = nvariables(equations)
    else
        data = convert_variables(u, solution_variables, mesh, equations, dg, cache)
        n_vars = size(data, 1)
    # else
    #     # Reinterpret the solution array as an array of conservative variables,
    #     # compute the solution variables via broadcasting, and reinterpret the
    #     # result as a plain array of floating point numbers
    #     data = Array(reinterpret(eltype(u),
    #                              solution_variables.(reinterpret(SVector{nvariables(equations),
    #                                                                      eltype(u)}, u),
    #                                                  Ref(equations))))

    #     # Find out variable count by looking at output from `solution_variables` function
    #     n_vars = size(data, 1)
    end

    # Open file (clobber existing content)
    # TODO: Create a function to do this in Trixi.jl to avoid duplicated code.
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