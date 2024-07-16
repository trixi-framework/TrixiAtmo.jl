@inline function reference_normal_vector(direction)
    # Get the normal vector to the reference element at a given node
    orientation = (direction + 1) >> 1
    if orientation == 1
        return (isodd(direction) ? -1 : 1) * SVector(1f0,0f0,0f0)
    else
        return (isodd(direction) ? -1 : 1) * SVector(0f0,1f0,0f0)
    end
end

@inline function cartesian2contravariant(u_node::SVector{NVARS}, 
                                         mesh::Trixi.AbstractMesh{2}, 
                                         equations::AbstractCovariantEquations3D{NVARS}, 
                                         i, j, element, cache) where {NVARS}

    (; contravariant_vectors, inverse_jacobian) = cache.elements

    Ja11, Ja12, Ja13 = Trixi.get_contravariant_vector(1, contravariant_vectors, i, j,
                                                      element)
    Ja21, Ja22, Ja23 = Trixi.get_contravariant_vector(2, contravariant_vectors, i, j, 
                                                      element)
    u_con_1, u_con_2 = inverse_jacobian[i,j,element] * 
            SMatrix{2,3}(Ja11, Ja21, Ja12, Ja22, Ja13, Ja23) * 
            SVector(u_node[1], u_node[2], u_node[3])

    return SVector{NVARS}(u_con_1, 
                          u_con_2, 
                          zero(eltype(u_node)),
                          u_node[4:end]...)
end

@inline function contravariant2cartesian(u_node::SVector{NVARS}, 
                                         mesh::Trixi.AbstractMesh{2}, 
                                         equations::AbstractCovariantEquations3D{NVARS}, 
                                         i, j, element, cache) where {NVARS}

    (; jacobian_matrix) = cache.elements

    u_car_1, u_car_2, u_car_3 = SMatrix{3,2}(jacobian_matrix[:,:, i,j,element]) * 
        SVector(u_node[1], u_node[2])

    return SVector{NVARS}(u_car_1, 
                          u_car_2, 
                          u_car_3, 
                          u_node[4:end]...)
end

function Trixi.compute_coefficients!(u, func, t, mesh::Trixi.AbstractMesh{2}, 
    equations::AbstractCovariantEquations3D{NVARS}, dg::DG, cache) where {NVARS}

    (; node_coordinates) = cache.elements

    for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            x_node = Trixi.get_node_coords(node_coordinates, equations, dg, i, j, element)
            
            u_car = func(x_node, t, equations)
            u_node = cartesian2contravariant(u_car, mesh, equations, i, j, element, cache)

            Trixi.set_node_vars!(u, u_node, equations, dg, i, j, element)
        end
    end
end

@inline function Trixi.weak_form_kernel!(du, u, element, 
    mesh::Union{StructuredMesh{2}, UnstructuredMesh2D, P4estMesh{2}, T8codeMesh{2}},
    nonconservative_terms::False, equations::AbstractCovariantEquations3D,
    dg::DGSEM, cache, alpha = true)
   
    (; derivative_dhat) = dg.basis
    (; inverse_jacobian) = cache.elements

    for j in eachnode(dg), i in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        
        contravariant_flux1 = flux(u_node, 1, equations) / inverse_jacobian[i,j,element]
        contravariant_flux2 = flux(u_node, 2, equations) / inverse_jacobian[i,j,element]

        for ii in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[ii, i],
                    contravariant_flux1, equations, dg, ii, j, element)
        end
        
        for jj in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[jj, j],
                    contravariant_flux2, equations, dg, i, jj, element)
        end
    end

    return nothing
end

function Trixi.calc_interface_flux!(surface_flux_values,
    mesh::Union{P4estMesh{2}, T8codeMesh{2}},
    nonconservative_terms,
    equations::AbstractCovariantEquations3D, surface_integral, dg::DG, cache)

    (; neighbor_ids, node_indices) = cache.interfaces
    index_range = eachnode(dg)
    index_end = last(index_range)

    for interface in Trixi.eachinterface(dg, cache)
        # Get element and side index information on the primary element
        primary_element = neighbor_ids[1, interface]
        primary_indices = node_indices[1, interface]
        primary_direction = Trixi.indices2direction(primary_indices)

        # Create the local i,j indexing
        i_primary_start, i_primary_step = Trixi.index_to_start_step_2d(primary_indices[1],
                                                                       index_range)
        j_primary_start, j_primary_step = Trixi.index_to_start_step_2d(primary_indices[2],
                                                                       index_range)

        i_primary = i_primary_start
        j_primary = j_primary_start

        # Get element and side index information on the secondary element
        secondary_element = neighbor_ids[2, interface]
        secondary_indices = node_indices[2, interface]
        secondary_direction = Trixi.indices2direction(secondary_indices)

        i_secondary_start, i_secondary_step = Trixi.index_to_start_step_2d(
            secondary_indices[1], index_range)
        j_secondary_start, j_secondary_step = Trixi.index_to_start_step_2d(
                secondary_indices[2], index_range)

        i_secondary = i_secondary_start
        j_secondary = j_secondary_start

        # Initiate the secondary index to be used in the surface for loop.
        # This index on the primary side will always run forward but
        # the secondary index might need to run backwards for flipped sides.
        if :i_backward in secondary_indices
            node_secondary = index_end
            node_secondary_step = -1
        else
            node_secondary = 1
            node_secondary_step = 1
        end

        for node in eachnode(dg)
                    
            Trixi.calc_interface_flux!(surface_flux_values, mesh, nonconservative_terms,
                                       equations, surface_integral, dg, cache, interface, 
                                       i_primary, j_primary, i_secondary, j_secondary,
                                       node, primary_direction, primary_element,
                                       node_secondary, secondary_direction,
                                       secondary_element)

            # Increment primary and secondary element indices
            i_primary += i_primary_step
            j_primary += j_primary_step

            i_secondary += i_secondary_step
            j_secondary += j_secondary_step

            # Increment the surface node index along the secondary element
            node_secondary += node_secondary_step
        end
    end

    return nothing
end

@inline function Trixi.calc_interface_flux!(surface_flux_values,
    mesh::Union{P4estMesh{2}, T8codeMesh{2}},
    nonconservative_terms::False, equations,
    surface_integral, dg::DG, cache, interface_index, 
    i_primary, j_primary, i_secondary, j_secondary,
    primary_node_index, primary_direction_index,
    primary_element_index,
    secondary_node_index, secondary_direction_index,
    secondary_element_index)

    (; inverse_jacobian) = cache.elements
    (; u) = cache.interfaces
    (; surface_flux) = surface_integral
    
    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, primary_node_index,
        interface_index)

    # transform to Cartesian components
    u_ll_car = contravariant2cartesian(u_ll, mesh, equations, i_primary, j_primary,
                                       primary_element_index, cache)
    u_rr_car = contravariant2cartesian(u_rr, mesh, equations, i_secondary, j_secondary,
                                       secondary_element_index, cache)
                                    
    # transform solution on secondary element into frame of primary element
    u_rr_ll = cartesian2contravariant(u_rr_car, mesh, equations, 
                                   i_primary, j_primary, primary_element_index, cache)

    # compute the flux for the primary element scaled by the Jacobian
    reference_normal_primary = reference_normal_vector(primary_direction_index)
    flux_primary = surface_flux(u_ll, u_rr_ll, reference_normal_primary, equations) / 
        inverse_jacobian[i_primary, j_primary, primary_element_index]

    # transform solution on primary element into frame of secondary element
    u_ll_rr = cartesian2contravariant(u_ll_car, mesh, equations, 
        i_secondary, j_secondary, secondary_element_index, cache)

    # compute the flux for the secondary element scaled by the Jacobian
    reference_normal_secondary = reference_normal_vector(secondary_direction_index)
    flux_secondary = surface_flux(u_rr, u_ll_rr, reference_normal_secondary, equations) / 
            inverse_jacobian[i_secondary, j_secondary, secondary_element_index]
    
    for v in eachvariable(equations)
        surface_flux_values[v, primary_node_index, primary_direction_index,
            primary_element_index] = flux_primary[v]
        surface_flux_values[v, secondary_node_index, secondary_direction_index,
            secondary_element_index] = flux_secondary[v]
    end
end
