@muladd begin
#! format: noindent

# Convert to another set of variables using a solution_variables function
function convert_variables(u, solution_variables, mesh::P4estMesh{2},
                           equations::AbstractEquations{3},
                           dg, cache)
    (; contravariant_vectors) = cache.elements
    # Extract the number of solution variables to be output 
    # (may be different than the number of conservative variables) 
    n_vars = length(Trixi.varnames(solution_variables, equations))

    # Allocate storage for output
    data = Array{eltype(u)}(undef, n_vars, nnodes(dg), nnodes(dg), nelements(dg, cache))

    # Loop over all nodes and convert variables, passing in auxiliary variables
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            normal_vector_node = Trixi.get_contravariant_vector(3,
                                                                contravariant_vectors,
                                                                i, j, element)

            if applicable(solution_variables, u, normal_vector_node, mesh, equations,
                          dg,
                          cache, i, j, element)
                # The solution_variables function depends on the solution on other nodes, the normal_vector_node, etc.
                data_node = solution_variables(u, normal_vector_node, mesh, equations,
                                               dg,
                                               cache, i, j, element)
            else
                # The solution_variables function depends on u_node and equations
                u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
                data_node = solution_variables(u_node, equations)
            end

            for v in 1:n_vars
                data[v, i, j, element] = data_node[v]
            end
        end
    end
    return data
end

function Trixi.save_solution_file(u, time, dt, timestep,
                                  mesh::P4estMesh{2},
                                  equations::Union{AbstractEquations{3},
                                                   AbstractCovariantEquations{2}},
                                  dg::DG, cache,
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

# Calculate the primitive variables and the relative vorticity at a given node
@inline function cons2prim_and_vorticity(u, normal_vector,
                                         mesh::P4estMesh{2},
                                         equations::AbstractEquations{3},
                                         dg::DGSEM, cache, i, j, element)

    # compute the vorticity and project onto normal vector
    vorticity = calc_vorticity_node(u, mesh, equations, dg, cache, i, j, element)
    relative_vorticity = dot(vorticity, normal_vector) / norm(normal_vector)

    u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)

    # Return the solution variables
    return SVector(cons2prim(u_node, equations)..., relative_vorticity)
end

@inline function calc_vorticity_node(u, mesh::P4estMesh{2},
                                     equations::AbstractEquations{3}, dg::DGSEM, cache,
                                     i, j, element)
    (; derivative_matrix) = dg.basis
    (; contravariant_vectors, inverse_jacobian) = cache.elements

    # Compute gradients in reference space
    dv1dxi1 = dv1dxi2 = zero(eltype(u))
    dv2dxi1 = dv2dxi2 = zero(eltype(u))
    dv3dxi1 = dv3dxi2 = zero(eltype(u))
    for ii in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, ii, j, element)
        v1, v2, v3 = velocity(u_node, equations)
        dv1dxi1 += derivative_matrix[i, ii] * v1
        dv2dxi1 += derivative_matrix[i, ii] * v2
        dv3dxi1 += derivative_matrix[i, ii] * v3
    end

    for jj in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, jj, element)
        v1, v2, v3 = velocity(u_node, equations)
        dv1dxi2 += derivative_matrix[j, jj] * v1
        dv2dxi2 += derivative_matrix[j, jj] * v2
        dv3dxi2 += derivative_matrix[j, jj] * v3
    end

    # Transform gradients to Cartesian space
    Ja11, Ja12, Ja13 = Trixi.get_contravariant_vector(1, contravariant_vectors, i, j,
                                                      element)
    Ja21, Ja22, Ja23 = Trixi.get_contravariant_vector(2, contravariant_vectors, i, j,
                                                      element)

    dv1dy = (Ja12 * dv1dxi1 + Ja22 * dv1dxi2) * inverse_jacobian[i, j, element]
    dv1dz = (Ja13 * dv1dxi1 + Ja23 * dv1dxi2) * inverse_jacobian[i, j, element]
    dv2dx = (Ja11 * dv2dxi1 + Ja21 * dv2dxi2) * inverse_jacobian[i, j, element]
    dv2dz = (Ja13 * dv2dxi1 + Ja23 * dv2dxi2) * inverse_jacobian[i, j, element]
    dv3dx = (Ja11 * dv3dxi1 + Ja21 * dv3dxi2) * inverse_jacobian[i, j, element]
    dv3dy = (Ja12 * dv3dxi1 + Ja22 * dv3dxi2) * inverse_jacobian[i, j, element]

    # compute the vorticity
    return SVector(dv3dy - dv2dz, dv1dz - dv3dx, dv2dx - dv1dy)
end

# Variable names for cons2prim_and_vorticity
function Trixi.varnames(::typeof(cons2prim_and_vorticity),
                        equations::ShallowWaterEquations3D)
    return (varnames(cons2prim, equations)..., "vorticity")
end
end # @muladd
