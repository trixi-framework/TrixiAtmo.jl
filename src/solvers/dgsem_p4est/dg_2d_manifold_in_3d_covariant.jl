@muladd begin
#! format: noindent

# Compute coefficients for an initial condition that uses auxiliary variables
function Trixi.compute_coefficients!(u, func, t, mesh::P4estMesh{2},
                                     equations::AbstractCovariantEquations{2}, dg::DG,
                                     cache)
    (; aux_node_vars) = cache.auxiliary_variables
    (; node_coordinates) = cache.elements

    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            # Get physical Cartesian coordinates at node (i, j)
            x_node = Trixi.get_node_coords(node_coordinates, equations, dg, i, j,
                                           element)

            # Get auxiliary variables at node (i, j)
            aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)

            # Compute nodal values of the provided function
            u_node = func(x_node, t, aux_node, equations)

            # Set nodal variables in cache
            Trixi.set_node_vars!(u, u_node, equations, dg, i, j, element)
        end
    end
end

# Weak form kernel which uses contravariant flux components, passing the geometric 
# information contained in the auxiliary variables to the flux function
@inline function Trixi.weak_form_kernel!(du, u, element, mesh::P4estMesh{2},
                                         nonconservative_terms::False,
                                         equations::AbstractCovariantEquations{2},
                                         dg::DGSEM, cache, alpha = true)
    (; derivative_dhat) = dg.basis
    (; aux_node_vars) = cache.auxiliary_variables

    for j in eachnode(dg), i in eachnode(dg)
        # Get solution variables and auxiliary variables at node (i, j)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)

        # Evaluate the contravariant flux components
        flux1 = flux(u_node, aux_node, 1, equations)
        flux2 = flux(u_node, aux_node, 2, equations)

        # Apply weak form derivative with respect to ξ¹ 
        for ii in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[ii, i],
                                             flux1, equations, dg, ii, j, element)
        end

        # Apply weak form derivative with respect to ξ²
        for jj in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[jj, j],
                                             flux2, equations, dg, i, jj, element)
        end
    end

    return nothing
end

# Flux differencing kernel which uses contravariant flux components, passing the geometric 
# information contained in the auxiliary variables to the flux function
@inline function Trixi.flux_differencing_kernel!(du, u, element, mesh::P4estMesh{2},
                                                 nonconservative_terms::False,
                                                 equations::AbstractCovariantEquations{2},
                                                 volume_flux, dg::DGSEM, cache,
                                                 alpha = true)
    (; derivative_split) = dg.basis
    (; aux_node_vars) = cache.auxiliary_variables

    for j in eachnode(dg), i in eachnode(dg)
        # Get solution variables and auxiliary variables at node (i, j)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)

        # ξ¹ direction
        for ii in (i + 1):nnodes(dg)
            # Get solution variables and auxiliary variables at node (ii, j)
            u_node_ii = Trixi.get_node_vars(u, equations, dg, ii, j, element)
            aux_node_ii = get_node_aux_vars(aux_node_vars, equations, dg, ii, j,
                                            element)

            # Evaluate contravariant component 1 of the variable-coefficient two-point flux
            flux1 = volume_flux(u_node, u_node_ii, aux_node, aux_node_ii, 1, equations)

            # Multiply by entry of split derivative matrix and add to right-hand side
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[i, ii], flux1,
                                             equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[ii, i], flux1,
                                             equations, dg, ii, j, element)
        end

        # ξ² direction
        for jj in (j + 1):nnodes(dg)
            # Get solution variables and auxiliary variables at node (i, jj)
            u_node_jj = Trixi.get_node_vars(u, equations, dg, i, jj, element)
            aux_node_jj = get_node_aux_vars(aux_node_vars, equations, dg, i, jj,
                                            element)

            # Evaluate contravariant component 2 of the variable-coefficient two-point flux
            flux2 = volume_flux(u_node, u_node_jj, aux_node, aux_node_jj, 2, equations)

            # Multiply by entry of split derivative matrix and add to right-hand side
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[j, jj], flux2,
                                             equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[jj, j], flux2,
                                             equations, dg, i, jj, element)
        end
    end

    return nothing
end

# Non-conservative flux differencing kernel which uses contravariant flux components, 
# passing the geometric information contained in the auxiliary variables to the flux 
# function
@inline function Trixi.flux_differencing_kernel!(du, u, element, mesh::P4estMesh{2},
                                                 nonconservative_terms::True,
                                                 equations::AbstractCovariantEquations{2},
                                                 volume_flux, dg::DGSEM, cache,
                                                 alpha = true)
    (; derivative_split) = dg.basis
    (; aux_node_vars) = cache.auxiliary_variables
    symmetric_flux, nonconservative_flux = volume_flux

    # Apply the symmetric flux as usual
    Trixi.flux_differencing_kernel!(du, u, element, mesh, False(), equations,
                                    symmetric_flux, dg, cache, alpha)

    for j in eachnode(dg), i in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)

        # ξ¹ direction
        integral_contribution = zero(u_node)
        for ii in eachnode(dg)
            u_node_ii = Trixi.get_node_vars(u, equations, dg, ii, j, element)
            aux_node_ii = get_node_aux_vars(aux_node_vars, equations, dg, ii, j,
                                            element)
            flux1 = nonconservative_flux(u_node, u_node_ii, aux_node, aux_node_ii, 1,
                                         equations)
            integral_contribution = integral_contribution +
                                    derivative_split[i, ii] * flux1
        end

        # ξ² direction
        for jj in eachnode(dg)
            u_node_jj = Trixi.get_node_vars(u, equations, dg, i, jj, element)
            aux_node_jj = get_node_aux_vars(aux_node_vars, equations, dg, i, jj,
                                            element)
            flux2 = nonconservative_flux(u_node, u_node_jj, aux_node, aux_node_jj, 2,
                                         equations)
            integral_contribution = integral_contribution +
                                    derivative_split[j, jj] * flux2
        end

        Trixi.multiply_add_to_node_vars!(du, alpha * 0.5f0, integral_contribution,
                                         equations, dg, i, j, element)
    end

    return nothing
end

# Pointwise interface flux in local coordinates for problems without nonconservative terms
@inline function Trixi.calc_interface_flux!(surface_flux_values, mesh::P4estMesh{2},
                                            nonconservative_terms::False,
                                            equations::AbstractCovariantEquations{2},
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

    # Get surface values for solution and auxiliary variables
    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, primary_node_index,
                                             interface_index)
    aux_vars_ll, aux_vars_rr = get_surface_node_aux_vars(aux_surface_node_vars,
                                                         equations,
                                                         dg, primary_node_index,
                                                         interface_index)

    # Compute flux in the primary element's coordinate system
    u_rr_global = contravariant2global(u_rr, aux_vars_rr, equations)
    u_rr_transformed_to_ll = global2contravariant(u_rr_global, aux_vars_ll,
                                                  equations)
    if isodd(primary_direction_index)
        flux_primary = -surface_flux(u_rr_transformed_to_ll, u_ll,
                                     aux_vars_ll, aux_vars_ll,
                                     (primary_direction_index + 1) >> 1, equations)
    else
        flux_primary = surface_flux(u_ll, u_rr_transformed_to_ll,
                                    aux_vars_ll, aux_vars_ll,
                                    (primary_direction_index + 1) >> 1, equations)
    end

    # Compute flux in the secondary element's coordinate system
    u_ll_global = contravariant2global(u_ll, aux_vars_ll, equations)
    u_ll_transformed_to_rr = global2contravariant(u_ll_global, aux_vars_rr,
                                                  equations)
    if isodd(secondary_direction_index)
        flux_secondary = -surface_flux(u_ll_transformed_to_rr, u_rr,
                                       aux_vars_rr, aux_vars_rr,
                                       (secondary_direction_index + 1) >> 1, equations)
    else
        flux_secondary = surface_flux(u_rr, u_ll_transformed_to_rr,
                                      aux_vars_rr, aux_vars_rr,
                                      (secondary_direction_index + 1) >> 1, equations)
    end

    # Update the surface flux values on both sides of the interface
    for v in eachvariable(equations)
        surface_flux_values[v, primary_node_index, primary_direction_index,
        primary_element_index] = flux_primary[v]
        surface_flux_values[v, secondary_node_index, secondary_direction_index,
        secondary_element_index] = flux_secondary[v]
    end
end

# Pointwise interface flux in local coordinates for problems with nonconservative terms
@inline function Trixi.calc_interface_flux!(surface_flux_values, mesh::P4estMesh{2},
                                            nonconservative_terms::True,
                                            equations::AbstractCovariantEquations{2},
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
    surface_flux, nonconservative_flux = surface_integral.surface_flux

    # Get surface values for solution and auxiliary variables
    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, primary_node_index,
                                             interface_index)
    aux_vars_ll, aux_vars_rr = get_surface_node_aux_vars(aux_surface_node_vars,
                                                         equations,
                                                         dg, primary_node_index,
                                                         interface_index)

    # Compute flux in the primary element's coordinate system
    u_rr_global = contravariant2global(u_rr, aux_vars_rr, equations)
    u_rr_transformed_to_ll = global2contravariant(u_rr_global, aux_vars_ll,
                                                  equations)
    primary_orientation = (primary_direction_index + 1) >> 1
    if isodd(primary_direction_index)
        flux_primary = -(surface_flux(u_rr_transformed_to_ll, u_ll, aux_vars_ll,
                                      aux_vars_ll, primary_orientation, equations) +
                         0.5f0 * nonconservative_flux(u_rr_transformed_to_ll, u_ll,
                                              aux_vars_ll, aux_vars_ll,
                                              primary_orientation, equations))
    else
        flux_primary = surface_flux(u_ll, u_rr_transformed_to_ll, aux_vars_ll,
                                    aux_vars_ll, primary_orientation, equations) +
                       0.5f0 * nonconservative_flux(u_ll, u_rr_transformed_to_ll,
                                            aux_vars_ll, aux_vars_ll,
                                            primary_orientation, equations)
    end

    # Compute flux in the secondary element's coordinate system
    u_ll_global = contravariant2global(u_ll, aux_vars_ll, equations)
    u_ll_transformed_to_rr = global2contravariant(u_ll_global, aux_vars_rr,
                                                  equations)
    secondary_orientation = (secondary_direction_index + 1) >> 1
    if isodd(secondary_direction_index)
        flux_secondary = -(surface_flux(u_ll_transformed_to_rr, u_rr, aux_vars_rr,
                                        aux_vars_rr, secondary_orientation, equations) +
                           0.5f0 * nonconservative_flux(u_ll_transformed_to_rr, u_rr,
                                                aux_vars_rr, aux_vars_rr,
                                                secondary_orientation, equations))
    else
        flux_secondary = surface_flux(u_rr, u_ll_transformed_to_rr, aux_vars_rr,
                                      aux_vars_rr, secondary_orientation, equations) +
                         0.5f0 * nonconservative_flux(u_rr, u_ll_transformed_to_rr,
                                              aux_vars_rr, aux_vars_rr,
                                              secondary_orientation, equations)
    end

    # Update the surface flux values on both sides of the interface
    for v in eachvariable(equations)
        surface_flux_values[v, primary_node_index, primary_direction_index,
        primary_element_index] = flux_primary[v]
        surface_flux_values[v, secondary_node_index, secondary_direction_index,
        secondary_element_index] = flux_secondary[v]
    end
end

# This function passes the auxiliary variables into the source term
function Trixi.calc_sources!(du, u, t, source_terms,
                             equations::AbstractCovariantEquations{2}, dg::DG, cache)
    (; aux_node_vars) = cache.auxiliary_variables

    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            u_local = Trixi.get_node_vars(u, equations, dg, i, j, element)
            x_local = Trixi.get_node_coords(cache.elements.node_coordinates, equations,
                                            dg,
                                            i, j, element)
            aux_local = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            du_local = source_terms(u_local, x_local, t, aux_local, equations)
            Trixi.add_to_node_vars!(du, du_local, equations, dg, i, j, element)
        end
    end

    return nothing
end

# Version of calc_sources! specialized for covariant formulation with no source term
function Trixi.calc_sources!(du, u, t, source_terms::Nothing,
                             equations::AbstractCovariantEquations{2}, dg::DG, cache)
    return nothing
end

# Apply the exact Jacobian stored in auxiliary variables
function Trixi.apply_jacobian!(du, mesh::P4estMesh{2},
                               equations::AbstractCovariantEquations{2},
                               dg::DG, cache)
    (; aux_node_vars) = cache.auxiliary_variables

    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            aux_node = get_node_aux_vars(aux_node_vars, equations, dg, i, j, element)
            factor = -1 / area_element(aux_node, equations)

            for v in eachvariable(equations)
                du[v, i, j, element] *= factor
            end
        end
    end

    return nothing
end
end # @muladd
