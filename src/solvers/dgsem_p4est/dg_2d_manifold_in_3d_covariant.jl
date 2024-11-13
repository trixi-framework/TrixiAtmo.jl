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

# Evaluate the initial condition in spherical vector components, then 
# transform to contravariant components for use as prognostic variables
function Trixi.compute_coefficients!(u, func, t, mesh::P4estMesh{2},
                                     equations::AbstractCovariantEquations{2}, dg::DG,
                                     cache)
    for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            x_node = Trixi.get_node_coords(cache.elements.node_coordinates,
                                           equations, dg, i, j, element)
            u_node = func(x_node, t, equations, cache, (i, j), element)
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

        contravariant_flux1 = flux(u_node, 1, equations, cache, (i, j), element)
        contravariant_flux2 = flux(u_node, 2, equations, cache, (i, j), element)

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

        # x direction
        for ii in (i + 1):nnodes(dg)
            u_node_ii = Trixi.get_node_vars(u, equations, dg, ii, j, element)
            flux1 = volume_flux(u_node, u_node_ii, 1, equations, cache,
                                (i, j), (ii, j), element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[i, ii], flux1,
                                             equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[ii, i], flux1,
                                             equations, dg, ii, j, element)
        end

        # y direction
        for jj in (j + 1):nnodes(dg)
            u_node_jj = Trixi.get_node_vars(u, equations, dg, i, jj, element)
            flux2 = volume_flux(u_node, u_node_jj, 2, equations, cache,
                                (i, j), (i, jj), element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[j, jj], flux2,
                                             equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[jj, j], flux2,
                                             equations, dg, i, jj, element)
        end
    end

    return nothing
end

# Interface flux which transforms the contravariant prognostic variables into the same 
# reference coordinate system on each side of the interface, then applies the numerical 
# flux in reference space, taking in the reference normal vector
function Trixi.calc_interface_flux!(surface_flux_values,
                                    mesh::P4estMesh{2},
                                    nonconservative_terms,
                                    equations::AbstractCovariantEquations{2},
                                    surface_integral, dg::DG, cache)
    (; neighbor_ids, node_indices) = cache.interfaces

    index_range = eachnode(dg)
    index_end = last(index_range)

    for interface in Trixi.eachinterface(dg, cache)

        # Get element and side index information on the primary element
        primary_element = neighbor_ids[1, interface]
        primary_indices = node_indices[1, interface]
        primary_direction = Trixi.indices2direction(primary_indices)

        # Create the local i,j indexing on the primary element
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

        # Create the local i,j indexing on the secondary element
        i_secondary_start, i_secondary_step = Trixi.index_to_start_step_2d(secondary_indices[1],
                                                                           index_range)
        j_secondary_start, j_secondary_step = Trixi.index_to_start_step_2d(secondary_indices[2],
                                                                           index_range)
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
                                       equations, surface_integral, dg, cache,
                                       interface, i_primary, j_primary, i_secondary,
                                       j_secondary,
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

# Pointwise interface flux, transforming the contravariant prognostic variables into the 
# local coordinate system
@inline function Trixi.calc_interface_flux!(surface_flux_values, mesh::P4estMesh{2},
                                            nonconservative_terms::False, equations,
                                            surface_integral, dg::DG, cache,
                                            interface_index,
                                            i_primary, j_primary,
                                            i_secondary, j_secondary,
                                            primary_node_index,
                                            primary_direction_index,
                                            primary_element_index,
                                            secondary_node_index,
                                            secondary_direction_index,
                                            secondary_element_index)
    (; u) = cache.interfaces
    (; surface_flux) = surface_integral

    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, primary_node_index,
                                             interface_index)
    # Gather multi-indices
    node_ll = (i_primary, j_primary)
    node_rr = (i_secondary, j_secondary)

    # Convert to spherical components on each element
    u_ll_spherical = contravariant2spherical(u_ll, equations, cache, node_ll,
                                             primary_element_index)
    u_rr_spherical = contravariant2spherical(u_rr, equations, cache,
                                             node_rr, secondary_element_index)

    # Evaluate u_rr in primary coordinate system 
    u_rr_transformed_to_ll = spherical2contravariant(u_rr_spherical, equations, cache,
                                                     node_ll, primary_element_index)
    # Compute flux on primary element
    if isodd(primary_direction_index)
        flux_primary = -surface_flux(u_rr_transformed_to_ll, u_ll,
                                     (primary_direction_index + 1) >> 1, equations,
                                     cache, node_ll, node_ll, primary_element_index)
    else
        flux_primary = surface_flux(u_ll, u_rr_transformed_to_ll,
                                    (primary_direction_index + 1) >> 1, equations,
                                    cache, node_ll, node_ll, primary_element_index)
    end
    # Evaluate u_ll in secondary coordinate system
    u_ll_transformed_to_rr = spherical2contravariant(u_ll_spherical, equations, cache,
                                                     node_rr, secondary_element_index)
    # Compute flux on secondary element
    if isodd(secondary_direction_index)
        flux_secondary = -surface_flux(u_ll_transformed_to_rr, u_rr,
                                       (secondary_direction_index + 1) >> 1, equations,
                                       cache, node_rr, node_rr, secondary_element_index)
    else
        flux_secondary = surface_flux(u_rr, u_ll_transformed_to_rr,
                                      (secondary_direction_index + 1) >> 1, equations,
                                      cache, node_rr, node_rr, secondary_element_index)
    end

    # Now we can update the surface flux values, where all vector variables are stored 
    # as contravariant components
    for v in eachvariable(equations)
        surface_flux_values[v, primary_node_index, primary_direction_index,
        primary_element_index] = flux_primary[v]
        surface_flux_values[v, secondary_node_index, secondary_direction_index,
        secondary_element_index] = flux_secondary[v]
    end
end
end # @muladd
