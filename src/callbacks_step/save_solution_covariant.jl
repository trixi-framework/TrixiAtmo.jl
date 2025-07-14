@muladd begin
#! format: noindent

# Convert to another set of variables using a solution_variables function
function convert_variables(u, solution_variables, mesh::P4estMesh{2},
                           equations::AbstractCovariantEquations{2}, dg, cache)
    (; aux_node_vars) = cache.auxiliary_variables
    # Extract the number of solution variables to be output 
    # (may be different than the number of conservative variables) 
    n_vars = length(Trixi.varnames(solution_variables, equations))

    # Allocate storage for output
    data = Array{eltype(u)}(undef, n_vars, nnodes(dg), nnodes(dg), nelements(dg, cache))

    # Loop over all nodes and convert variables, passing in auxiliary variables
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            if applicable(solution_variables, u, equations, dg, cache, i, j, element)
                # The solution_variables function depends on the solution on other nodes
                data_node = solution_variables(u, equations, dg, cache, i, j, element)
            else
                u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
                aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j,
                                             element)
                data_node = solution_variables(u_node, aux_node, equations)
            end

            for v in 1:n_vars
                data[v, i, j, element] = data_node[v]
            end
        end
    end
    return data
end

# Convert to another set of variables using a solution_variables function
function convert_variables(u, solution_variables, mesh::DGMultiMesh,
                           equations::AbstractCovariantEquations, dg, cache)
    (; aux_values) = cache
    # Extract the number of solution variables to be output 
    # (may be different than the number of conservative variables) 
    n_vars = length(Trixi.varnames(solution_variables, equations))

    # Allocate storage for output, an Np x n_elements array of nvars arrays
    # data = Array{eltype(u)}(undef, n_vars, (size(aux_node_vars)[2:end])...)
    data = map(_ -> SVector{n_vars, eltype(u)}(zeros(n_vars)), u)
    
    # Loop over all nodes and convert variables, passing in auxiliary variables
    for i in Trixi.each_dof_global(mesh, dg)
        u_node = u[:, i]
        aux_node = aux_values[i]
        data[i] = solution_variables(u_node, aux_node, equations)
    end
    return data
end

# Calculate the primitive variables and relative vorticity at a given node. The velocity
# components in the global coordinate system and the bottom topography are returned, such 
# that the outputs for CovariantShallowWaterEquations2D and ShallowWaterEquations3D are the 
# same variables, provided that GlobalCartesianCoordinates are used in the former case.
@inline function cons2prim_and_vorticity(u, equations::AbstractCovariantEquations{2},
                                         dg::DGSEM, cache, i, j, element)
    (; aux_node_vars) = cache.auxiliary_variables
    u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
    aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
    relative_vorticity = calc_vorticity_node(u, equations, dg, cache, i, j, element)
    h_s = bottom_topography(aux_node, equations)
    primitive_global = contravariant2global(cons2prim(u_node, aux_node, equations),
                                            aux_node, equations)
    return SVector(primitive_global..., h_s, relative_vorticity)
end

# Calculate relative vorticity ζ = (∂₁v₂ - ∂₂v₁)/J for equations in covariant form
@inline function calc_vorticity_node(u, equations::AbstractCovariantEquations{2},
                                     dg::DGSEM, cache, i, j, element)
    (; derivative_matrix) = dg.basis
    (; aux_node_vars) = cache.auxiliary_variables

    dv2dxi1 = dv1dxi2 = zero(eltype(u))
    for ii in eachnode(dg)
        u_node_ii = Trixi.get_node_vars(u, equations, dg, ii, j, element)
        aux_node_ii = get_node_aux_vars(aux_node_vars, equations, dg, ii, j, element)
        vcov = metric_covariant(aux_node_ii, equations) *
               velocity_contravariant(u_node_ii, equations)
        dv2dxi1 = dv2dxi1 + derivative_matrix[i, ii] * vcov[2]
    end

    for jj in eachnode(dg)
        u_node_jj = Trixi.get_node_vars(u, equations, dg, i, jj, element)
        aux_node_jj = get_node_aux_vars(aux_node_vars, equations, dg, i, jj, element)
        vcov = metric_covariant(aux_node_jj, equations) *
               velocity_contravariant(u_node_jj, equations)
        dv1dxi2 = dv1dxi2 + derivative_matrix[j, jj] * vcov[1]
    end

    # compute the relative vorticity
    aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
    return (dv2dxi1 - dv1dxi2) / area_element(aux_node, equations)
end

# Variable names for cons2prim_and_vorticity, chosen to match those from 
# ShallowWaterEquations3D
function Trixi.varnames(::typeof(cons2prim_and_vorticity),
                        equations::AbstractCovariantShallowWaterEquations2D)
    return ("H", "v1", "v2", "v3", "b", "vorticity")
end
end # @muladd
