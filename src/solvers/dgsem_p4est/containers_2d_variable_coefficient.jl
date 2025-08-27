@muladd begin

# Container for storing values of auxiliary variables at volume/surface quadrature nodes
struct P4estAuxNodeVarsContainer{NDIMS, uEltype <: Real, NDIMSP2, NDIMSP3, AuxField}
    aux_node_vars::Array{uEltype, NDIMSP2} # [variable, i, j, k, element]
    aux_surface_node_vars::Array{uEltype, NDIMSP2} # [variable, i, j, k, element]
    aux_boundary_node_vars::Array{uEltype, NDIMSP2} # [leftright, var, i, boundary]
    aux_mortar_node_vars::Array{uEltype, NDIMSP3}   # [leftright, var, pos, i, mortar]


    # internal `resize!`able storage
    _aux_node_vars::Vector{uEltype}
    _aux_surface_node_vars::Vector{uEltype}
    _aux_boundary_node_vars::Vector{uEltype}
    _aux_mortar_node_vars::Vector{uEltype}

    # auxiliary field 
    aux_field::AuxField

end



function init_aux_vars(mesh::P4estMesh{2},
                        equations::AbstractVariableCoefficientEquations{2,1}, 
                        solver, cache, aux_field)
    @unpack elements, interfaces, boundaries, mortars = cache

    nelements = Trixi.ncells(mesh)
    ninterfaces = Trixi.count_required_surfaces(mesh).interfaces
    nboundaries = Trixi.nboundaries(boundaries)
    nmortars = Trixi.nmortars(mortars)
    NDIMS = Trixi.ndims(elements)
    uEltype = eltype(elements)


    _aux_node_vars = Vector{uEltype}(undef,
                                     n_aux_node_vars(equations) *
                                     nnodes(solver)^NDIMS * nelements)
    aux_node_vars = Trixi.unsafe_wrap(Array, pointer(_aux_node_vars),
                                      (n_aux_node_vars(equations),
                                       ntuple(_ -> nnodes(solver), NDIMS)...,
                                       nelements))

    _aux_surface_node_vars = Vector{uEltype}(undef,
                                             2 * n_aux_node_vars(equations) *
                                             nnodes(solver)^(NDIMS - 1) *
                                             ninterfaces)
    aux_surface_node_vars = Trixi.unsafe_wrap(Array,
                                              pointer(_aux_surface_node_vars),
                                              (2, n_aux_node_vars(equations),
                                               ntuple(_ -> nnodes(solver),
                                                      NDIMS - 1)...,
                                               ninterfaces))

    _aux_boundary_node_vars = Vector{uEltype}(undef,2 * n_aux_node_vars(equations) *
                                   nnodes(solver)^(NDIMS - 1) *
                                   nboundaries)

    aux_boundary_node_vars = unsafe_wrap(Array,
                                         pointer(_aux_boundary_node_vars),
                                         (2, n_aux_node_vars(equations),
                                          ntuple(_ -> nnodes(solver),
                                                 NDIMS - 1)...,
                                          nboundaries))

    _aux_mortar_node_vars = Vector{uEltype}(undef, 2 * n_aux_node_vars(equations) *
                                 2^(NDIMS - 1) * nnodes(solver)^(NDIMS - 1) *
                                 nmortars)

    aux_mortar_node_vars = unsafe_wrap(Array,
                                       pointer(_aux_mortar_node_vars),
                                       (2, n_aux_node_vars(equations),
                                        2^(NDIMS - 1),
                                        ntuple(_ -> nnodes(solver),
                                               NDIMS - 1)...,
                                        nmortars))

    aux_vars = P4estAuxNodeVarsContainer{NDIMS, uEltype, NDIMS + 2, NDIMS + 3, typeof(aux_field)}(aux_node_vars, aux_surface_node_vars,
                                        aux_boundary_node_vars, aux_mortar_node_vars,
                                        _aux_node_vars, _aux_surface_node_vars, 
                                        _aux_boundary_node_vars, _aux_mortar_node_vars, aux_field)

    init_aux_vars!(aux_vars, mesh, equations, solver, cache)

    return aux_vars
end

function init_aux_vars!(aux_vars, mesh, equations, solver, cache)
    init_aux_node_vars!(aux_vars, mesh, equations, solver, cache)
    init_aux_surface_node_vars!(aux_vars, mesh, equations, solver, cache)
    init_aux_boundary_node_vars!(aux_vars, mesh, equations, solver, cache)
    init_aux_mortar_node_vars!(aux_vars, mesh, equations, solver, cache)
end

#initialize the auxiliary variables at the nodes
function init_aux_node_vars!(
    aux_vars::P4estAuxNodeVarsContainer{2},
    mesh::P4estMesh{2},
    equations::AbstractVariableCoefficientEquations{2,1},
    solver,
    cache)

    (; node_coordinates) = cache.elements
    (; aux_node_vars, aux_field) = aux_vars


    Trixi.@threaded for element in eachelement(solver, cache) #go through the elements 
        for j in eachnode(solver), i in eachnode(solver) #and the nodes in both directions in the elements 
            if !isnothing(aux_field)
                x_node = Trixi.get_node_coords(node_coordinates, equations, solver, i, j, element)
                aux_node_vars[1, i, j, element] = aux_field(x_node)[1]
                aux_node_vars[2, i, j, element] = aux_field(x_node)[2]
            else
                aux_node_vars[1, i, j, element] = zero(eltype(aux_node_vars))
                aux_node_vars[2, i, j, element] = zero(eltype(aux_node_vars))
            end
        end
    end
    return nothing
end


# initialize auxiliary surface node variables
function init_aux_surface_node_vars!(
    aux_vars,
    mesh::P4estMesh{2},
    equations::AbstractVariableCoefficientEquations{2,1},
    solver,
    cache)

    (; neighbor_ids, node_indices) = cache.interfaces
    (; aux_node_vars, aux_surface_node_vars) = aux_vars
    
    index_range = eachnode(solver)


    Trixi.@threaded for interface in Trixi.eachinterface(solver, cache)
        primary_element = neighbor_ids[1, interface]
        secondary_element = neighbor_ids[2, interface]

        primary_indices = node_indices[1, interface]
        secondary_indices = node_indices[2, interface]
    
        i_p_start, i_p_step = Trixi.index_to_start_step_2d(primary_indices[1], index_range)
        j_p_start, j_p_step = Trixi.index_to_start_step_2d(primary_indices[2], index_range)
        
        i_s_start, i_s_step = Trixi.index_to_start_step_2d(secondary_indices[1], index_range)
        j_s_start, j_s_step = Trixi.index_to_start_step_2d(secondary_indices[2], index_range)
        
        #set startpoints to loop over the nodes for primary and secondary elements
        i_p = i_p_start
        j_p = j_p_start
        i_s = i_s_start
        j_s = j_s_start

        for node in eachnode(solver)
            for v in axes(aux_surface_node_vars, 2)
                aux_surface_node_vars[1, v, node, interface] =
                    aux_node_vars[v, i_p, j_p, primary_element]

                aux_surface_node_vars[2, v, node, interface] =
                    aux_node_vars[v, i_s, j_s, secondary_element]
            end

            # Step along face for both elements
            i_p += i_p_step
            j_p += j_p_step
            i_s += i_s_step
            j_s += j_s_step
        end 
    end
    return nothing
end

function init_aux_boundary_node_vars!(
    aux_vars,
    mesh::P4estMesh{2},
    equations::AbstractVariableCoefficientEquations{2,1},
    solver,
    cache)

    (; neighbor_ids, node_indices) = cache.boundaries
    (; aux_node_vars, aux_boundary_node_vars) = aux_vars
    
    index_range = eachnode(solver)

    Trixi.@threaded for boundary in Trixi.eachboundary(solver, cache)  # Loop over each boundary element 
        element = neighbor_ids[boundary]  # Get the element adjacent to this boundary
        
        face_indices = node_indices[boundary] 

        i_start, i_step = Trixi.index_to_start_step_2d(face_indices[1], index_range)
        j_start, j_step = Trixi.index_to_start_step_2d(face_indices[2], index_range)

        i = i_start
        j = j_start


        for l in eachnode(solver), v in axes(aux_boundary_node_vars, 2)
                    aux_boundary_node_vars[1, v, l, boundary] = aux_node_vars[v,
                                                                              i,
                                                                              j,
                                                                              element] 
        end

        i += i_step
        j += j_step
    end
    return nothing
end

function init_aux_mortar_node_vars!(
    aux_vars,
    mesh::P4estMesh{2},
    equations::AbstractVariableCoefficientEquations{2,1},
    solver,
    cache)

    (; neighbor_ids, node_indices) = cache.mortars
    (; aux_node_vars, aux_mortar_node_vars) = aux_vars

    index_range = eachnode(solver)


    Trixi.@threaded for mortar in Trixi.eachmortar(solver, cache) # Loop over each mortar
        # Get the IDs of both small neighbor elements on the mortar interface
        first_element = neighbor_ids[1, mortar]
        second_element = neighbor_ids[2, mortar]

        first_face_indices = node_indices[1, mortar]
        second_face_indices = node_indices[2, mortar]

        i1_start, i1_step = Trixi.index_to_start_step_2d(first_face_indices[1], index_range)
        j1_start, j1_step = Trixi.index_to_start_step_2d(first_face_indices[2], index_range)

        i2_start, i2_step = Trixi.index_to_start_step_2d(second_face_indices[1], index_range)
        j2_start, j2_step = Trixi.index_to_start_step_2d(second_face_indices[2], index_range)

        i1 = i1_start
        j1 = j1_start

        i2 = i2_start
        j2 = j2_start

        for l in eachnode(solver), v in axes(aux_mortar_node_vars, 2)
            aux_mortar_node_vars[:, v, 1, l, mortar] .= aux_node_vars[v, i1, j1, first_element]
            aux_mortar_node_vars[:, v, 2, l, mortar] .= aux_node_vars[v, i2, j2, second_element]
        
            i1 += i1_step
            j1 += j1_step

            i2 += i2_step
            j2 += j2_step
        end
    end
    return nothing
end

function Trixi.init_elements(
    mesh::P4estMesh{2,2,RealT},
    equations::VariableCoefficientAdvectionEquation2D,
    basis,
    metric_terms,
    ::Type{uEltype},
) where {RealT<:Real,uEltype<:Real}
    nelements = Trixi.ncells(mesh)

    _node_coordinates = Vector{RealT}(undef, 2 * nnodes(basis)^2 * nelements)
    node_coordinates = unsafe_wrap(
        Array,
        pointer(_node_coordinates),
        (2, ntuple(_ -> nnodes(basis), 2)..., nelements),
    )

    _jacobian_matrix = Vector{RealT}(undef, 2^2 * nnodes(basis)^2 * nelements)
    jacobian_matrix = unsafe_wrap(
        Array,
        pointer(_jacobian_matrix),
        (2, 2, ntuple(_ -> nnodes(basis), 2)..., nelements),
    )

    _contravariant_vectors = similar(_jacobian_matrix)
    contravariant_vectors =
        unsafe_wrap(Array, pointer(_contravariant_vectors), size(jacobian_matrix))

    _inverse_jacobian = Vector{RealT}(undef, nnodes(basis)^2 * nelements)
    inverse_jacobian = unsafe_wrap(
        Array,
        pointer(_inverse_jacobian),
        (ntuple(_ -> nnodes(basis), 2)..., nelements),
    )

    _surface_flux_values = Vector{uEltype}(
        undef,
        nvariables(equations) * nnodes(basis)^(2 - 1) * (2 * 2) * nelements,
    )
    surface_flux_values = unsafe_wrap(
        Array,
        pointer(_surface_flux_values),
        (
            nvariables(equations),
            ntuple(_ -> nnodes(basis), 2 - 1)...,
            2 * 2,
            nelements,
        ),
    )

    elements = Trixi.P4estElementContainer{2,RealT,uEltype,2 + 1,2 + 2,2 + 3}(
        node_coordinates,
        jacobian_matrix,
        contravariant_vectors,
        inverse_jacobian,
        surface_flux_values,
        _node_coordinates,
        _jacobian_matrix,
        _contravariant_vectors,
        _inverse_jacobian,
        _surface_flux_values,
    )

    Trixi.init_elements!(elements, mesh, basis)
    return elements

end #muladd
end
