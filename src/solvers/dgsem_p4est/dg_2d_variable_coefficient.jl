@inline function Trixi.weak_form_kernel!(du, u, 
                                   element, 
                                   mesh::P4estMesh{2},
                                   nonconservative_terms::Trixi.False,
                                   equation::AbstractVariableCoefficientEquations{2},
                                   dg::DGSEM, cache, alpha = true)  
    (; derivative_dhat) = dg.basis
    (; contravariant_vectors) = cache.elements
    (; aux_node_vars) = cache.auxiliary_variables

    for j in eachnode(dg), i in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)


        flux1 = flux(u_node, aux_node, 1, equations)
        flux2 = flux(u_node, aux_node, 2, equations)

        # Compute the contravariant flux by taking the scalar product of the
        # first contravariant vector Ja^1 and the flux vector
        Ja11, Ja12 = Trixi.get_contravariant_vector(1, contravariant_vectors, i, j, element)
        contravariant_flux1 = Ja11 * flux1 + Ja12 * flux2
        for ii in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[ii, i],
                                       contravariant_flux1, equations, dg, ii, j,
                                       element)
        end

        # Compute the contravariant flux by taking the scalar product of the
        # second contravariant vector Ja^2 and the flux vector
        Ja21, Ja22 = Trixi.get_contravariant_vector(2, contravariant_vectors, i, j, element)
        contravariant_flux2 = Ja21 * flux1 + Ja22 * flux2
        for jj in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[jj, j],
                                       contravariant_flux2, equations, dg, i, jj,
                                       element)
        end
    end
end



@inline function Trixi.calc_interface_flux!(surface_flux_values, mesh::P4estMesh{2},
    nonconservative_terms::Trixi.False,
    equations::AbstractVariableCoefficientEquations{2},
    surface_integral, dg::DG, cache,
    interface_index, normal_direction,
    primary_node_index,
    primary_direction_index,
    primary_element_index,
    secondary_node_index,
    secondary_direction_index,
    secondary_element_index)

    (; u) = cache.interfaces
    (; aux_surface_node_vars) = cache.auxiliary_variables
    (; surface_flux) = surface_integral

    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, primary_node_index,
                                       interface_index)

    # Get auxiliary variables
    aux_vars_ll, aux_vars_rr = get_surface_node_aux_vars(aux_surface_node_vars,
                                                         equations,
                                                         dg, primary_node_index,
                                                         interface_index)
    
                                                         #compute flux 
    flux_ = surface_flux(u_ll, u_rr, normal_direction, equations, aux_vars_ll, aux_vars_rr)
    
    # Store flux values
    for v in eachvariable(equations)
        surface_flux_values[v, primary_node_index, primary_direction_index, primary_element_index] = flux_[v]
        surface_flux_values[v, secondary_node_index, secondary_direction_index, secondary_element_index] = -flux_[v]
    end
end 
