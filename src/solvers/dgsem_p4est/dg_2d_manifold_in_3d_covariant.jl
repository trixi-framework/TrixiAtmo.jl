@muladd begin
#! format: noindent

# Custom rhs! for the covariant form. Note that the mortar flux is not included, as we do
# not yet support non-conforming meshes for the covariant solver.
function Trixi.rhs!(du, u, t, mesh::P4estMesh{2},
                    equations::AbstractCovariantEquations{2},
                    boundary_conditions, source_terms::Source,
                    dg::DG, cache) where {Source}
    # Reset du
    Trixi.@trixi_timeit Trixi.timer() "reset ∂u/∂t" Trixi.reset_du!(du, dg, cache)

    # Calculate volume integral
    Trixi.@trixi_timeit Trixi.timer() "volume integral" begin
        Trixi.calc_volume_integral!(du, u, mesh,
                                    Trixi.have_nonconservative_terms(equations),
                                    equations,
                                    dg.volume_integral, dg, cache)
    end

    # Prolong solution to interfaces
    Trixi.@trixi_timeit Trixi.timer() "prolong2interfaces" begin
        Trixi.prolong2interfaces!(cache, u, mesh, equations,
                                  dg.surface_integral, dg)
    end

    # Calculate interface fluxes
    Trixi.@trixi_timeit Trixi.timer() "interface flux" begin
        Trixi.calc_interface_flux!(cache.elements.surface_flux_values, mesh,
                                   Trixi.have_nonconservative_terms(equations),
                                   equations,
                                   dg.surface_integral, dg, cache)
    end

    # Prolong solution to boundaries
    Trixi.@trixi_timeit Trixi.timer() "prolong2boundaries" begin
        Trixi.prolong2boundaries!(cache, u, mesh, equations,
                                  dg.surface_integral, dg)
    end

    # Calculate boundary fluxes
    Trixi.@trixi_timeit Trixi.timer() "boundary flux" begin
        Trixi.calc_boundary_flux!(cache, t, boundary_conditions, mesh, equations,
                                  dg.surface_integral, dg)
    end

    # Calculate surface integrals
    Trixi.@trixi_timeit Trixi.timer() "surface integral" begin
        Trixi.calc_surface_integral!(du, u, mesh, equations,
                                     dg.surface_integral, dg, cache)
    end

    # Apply Jacobian from mapping to reference element
    Trixi.@trixi_timeit Trixi.timer() "Jacobian" Trixi.apply_jacobian!(du, mesh,
                                                                       equations, dg,
                                                                       cache)

    # Calculate source terms
    Trixi.@trixi_timeit Trixi.timer() "source terms" begin
        Trixi.calc_sources!(du, u, t, source_terms, equations, dg, cache)
    end

    return nothing
end

# Evaluate the initial condition in spherical vector components, then transform to 
# contravariant components for use as prognostic variables.
function Trixi.compute_coefficients!(u, func, t, mesh::P4estMesh{2},
                                     equations::AbstractCovariantEquations{2}, dg::DG,
                                     cache)

    # Store the auxiliary variables at the volume quadrature nodes in cache.elements
    compute_auxiliary_variables!(cache.elements, mesh, equations, dg)

    # Copy the appropriate interface values of the auxiliary variables into cache.interfaces
    prolong2interfaces_auxiliary!(cache, mesh, dg)

    for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            x_node = Trixi.get_node_coords(cache.elements.node_coordinates,
                                           equations, dg, i, j, element)
            aux_vars_node = get_node_aux_vars(cache.elements.auxiliary_variables,
                                              equations, dg, i, j, element)
            u_node = func(x_node, t, aux_vars_node, equations)
            Trixi.set_node_vars!(u, u_node, equations, dg, i, j, element)
        end
    end
end

# Weak form kernel which uses contravariant flux components, passing the element container 
# and node/element indices into the flux function to give the flux access to geometric 
# quantities
@inline function Trixi.weak_form_kernel!(du, u, element, mesh::P4estMesh{2},
                                         nonconservative_terms::False,
                                         equations::AbstractCovariantEquations{2},
                                         dg::DGSEM, cache, alpha = true)
    (; derivative_dhat) = dg.basis

    for j in eachnode(dg), i in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)

        aux_vars_node = get_node_aux_vars(cache.elements.auxiliary_variables,
                                          equations, dg, i, j, element)
        contravariant_flux1 = flux(u_node, aux_vars_node, 1, equations)
        contravariant_flux2 = flux(u_node, aux_vars_node, 2, equations)

        for ii in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[ii, i],
                                             contravariant_flux1, equations, dg, ii, j,
                                             element)
        end

        for jj in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_dhat[jj, j],
                                             contravariant_flux2, equations, dg, i, jj,
                                             element)
        end
    end

    return nothing
end

# Flux differencing kernel which uses contravariant flux components, passing the element
# container and node/element indices into the two-point volume flux function to give the 
# flux access to geometric quantities
@inline function Trixi.flux_differencing_kernel!(du, u, element, mesh::P4estMesh{2},
                                                 nonconservative_terms::False,
                                                 equations::AbstractCovariantEquations{2},
                                                 volume_flux, dg::DGSEM, cache,
                                                 alpha = true)
    (; derivative_split) = dg.basis

    for j in eachnode(dg), i in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        aux_vars_node = get_node_aux_vars(cache.elements.auxiliary_variables,
                                          equations, dg, i, j, element)

        # x direction
        for ii in (i + 1):nnodes(dg)
            u_node_ii = Trixi.get_node_vars(u, equations, dg, ii, j, element)
            aux_vars_node_ii = get_node_aux_vars(cache.elements.auxiliary_variables,
                                                 equations, dg, ii, j, element)

            flux1 = volume_flux(u_node, u_node_ii,
                                aux_vars_node, aux_vars_node_ii, 1, equations)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[i, ii], flux1,
                                             equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[ii, i], flux1,
                                             equations, dg, ii, j, element)
        end

        # y direction
        for jj in (j + 1):nnodes(dg)
            u_node_jj = Trixi.get_node_vars(u, equations, dg, i, jj, element)
            aux_vars_node_jj = get_node_aux_vars(cache.elements.auxiliary_variables,
                                                 equations, dg, i, jj, element)

            flux2 = volume_flux(u_node, u_node_jj,
                                aux_vars_node, aux_vars_node_jj, 2, equations)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[j, jj], flux2,
                                             equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[jj, j], flux2,
                                             equations, dg, i, jj, element)
        end
    end

    return nothing
end

# Pointwise interface flux, transforming the contravariant prognostic variables into the 
# local coordinate system
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
    (; u, auxiliary_variables) = cache.interfaces
    (; surface_flux) = surface_integral

    # Get surface values for solution and auxiliary variables
    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, primary_node_index,
                                             interface_index)
    aux_vars_ll, aux_vars_rr = get_surface_node_aux_vars(auxiliary_variables, equations,
                                                         dg, primary_node_index,
                                                         interface_index)

    # Compute flux in the primary element's coordinate system
    u_rr_spherical = contravariant2spherical(u_rr, aux_vars_rr, equations)
    u_rr_transformed_to_ll = spherical2contravariant(u_rr_spherical, aux_vars_ll,
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
    u_ll_spherical = contravariant2spherical(u_ll, aux_vars_ll, equations)
    u_ll_transformed_to_rr = spherical2contravariant(u_ll_spherical, aux_vars_rr,
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

# Apply the exact Jacobian stored in auxiliary variables
function Trixi.apply_jacobian!(du, mesh::P4estMesh{2},
                               equations::AbstractCovariantEquations{2},
                               dg::DG, cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            aux_vars_node = get_node_aux_vars(cache.elements.auxiliary_variables,
                                              equations, dg, i, j, element)
            factor = -1 / volume_element(aux_vars_node, equations)

            for v in eachvariable(equations)
                du[v, i, j, element] *= factor
            end
        end
    end

    return nothing
end

# Prolong the auxiliary variables and store in cache.interfaces
function prolong2interfaces_auxiliary!(cache, mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                       dg::DG)
    (; elements, interfaces) = cache
    index_range = eachnode(dg)

    Trixi.@threaded for interface in Trixi.eachinterface(dg, cache)
        # Copy solution data from the primary element using "delayed indexing" with
        # a start value and a step size to get the correct face and orientation.
        # Note that in the current implementation, the interface will be
        # "aligned at the primary element", i.e., the index of the primary side
        # will always run forwards.
        primary_element = interfaces.neighbor_ids[1, interface]
        primary_indices = interfaces.node_indices[1, interface]

        i_primary_start, i_primary_step = Trixi.index_to_start_step_2d(primary_indices[1],
                                                                       index_range)
        j_primary_start, j_primary_step = Trixi.index_to_start_step_2d(primary_indices[2],
                                                                       index_range)

        i_primary = i_primary_start
        j_primary = j_primary_start
        for i in eachnode(dg)
            for v in axes(elements.auxiliary_variables, 1)
                interfaces.auxiliary_variables[1, v, i, interface] = elements.auxiliary_variables[v,
                                                                                                  i_primary,
                                                                                                  j_primary,
                                                                                                  primary_element]
            end
            i_primary += i_primary_step
            j_primary += j_primary_step
        end

        # Copy solution data from the secondary element using "delayed indexing" with
        # a start value and a step size to get the correct face and orientation.
        secondary_element = interfaces.neighbor_ids[2, interface]
        secondary_indices = interfaces.node_indices[2, interface]

        i_secondary_start, i_secondary_step = Trixi.index_to_start_step_2d(secondary_indices[1],
                                                                           index_range)
        j_secondary_start, j_secondary_step = Trixi.index_to_start_step_2d(secondary_indices[2],
                                                                           index_range)

        i_secondary = i_secondary_start
        j_secondary = j_secondary_start
        for i in eachnode(dg)
            for v in axes(elements.auxiliary_variables, 1)
                interfaces.auxiliary_variables[2, v, i, interface] = elements.auxiliary_variables[v,
                                                                                                  i_secondary,
                                                                                                  j_secondary,
                                                                                                  secondary_element]
            end
            i_secondary += i_secondary_step
            j_secondary += j_secondary_step
        end
    end

    return nothing
end
end # @muladd
