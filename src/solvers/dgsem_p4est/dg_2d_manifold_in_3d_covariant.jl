@muladd begin
#! format: noindent

# Custom rhs! for the covariant form. Note that the mortar flux is not included, as we do
# not yet support non-conforming meshes for the covariant solver.
function Trixi.rhs!(du, u, t, mesh::P4estMesh{2},
                    equations::AbstractCovariantEquations{2},
                    initial_condition, boundary_conditions, source_terms::Source,
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
function Trixi.compute_coefficients!(u, func, t, mesh::P4estMesh,
                                     equations::AbstractCovariantEquations, dg::DG,
                                     cache)
    for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            x_node = Trixi.get_node_coords(cache.elements.node_coordinates,
                                           equations, dg, i, j, element)
            u_spherical = func(x_node, t, equations)
            u_con = spherical2contravariant(u_spherical, equations,
                                            cache.elements, i, j, element)
            Trixi.set_node_vars!(u, u_con, equations, dg, i, j, element)
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

        contravariant_flux1 = flux(u_node, 1, equations, cache.elements, i, j, element)
        contravariant_flux2 = flux(u_node, 2, equations, cache.elements, i, j, element)

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
            flux1 = volume_flux(u_node, u_node_ii, 1, equations,
                                cache.elements, i, j, ii, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[i, ii], flux1,
                                             equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[ii, i], flux1,
                                             equations, dg, ii, j, element)
        end

        # y direction
        for jj in (j + 1):nnodes(dg)
            u_node_jj = Trixi.get_node_vars(u, equations, dg, i, jj, element)
            flux2 = volume_flux(u_node, u_node_jj, 2, equations,
                                cache.elements, i, j, i, jj, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[j, jj], flux2,
                                             equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[jj, j], flux2,
                                             equations, dg, i, jj, element)
        end
    end

    return nothing
end

# Outward unit normal vector in reference coordinates for a quadrilateral element
@inline function reference_normal_vector(direction)
    orientation = (direction + 1) >> 1
    sign = isodd(direction) ? -1 : 1

    if orientation == 1
        return SVector(sign, 0.0f0)
    else
        return SVector(0.0f0, sign)
    end
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

    Trixi.@threaded for interface in Trixi.eachinterface(dg, cache)

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

    # Convert to spherical components on each element
    u_ll_pol = contravariant2spherical(u_ll, equations, cache.elements,
                                       i_primary, j_primary, primary_element_index)
    u_rr_pol = contravariant2spherical(u_rr, equations, cache.elements,
                                       i_secondary, j_secondary,
                                       secondary_element_index)

    # evaluate u_rr in secondary coordinate system 
    u_rr_ll = spherical2contravariant(u_rr_pol, equations, cache.elements,
                                      i_primary, j_primary, primary_element_index)

    # evaluate u_ll in primary coordinate system
    u_ll_rr = spherical2contravariant(u_ll_pol, equations, cache.elements,
                                      i_secondary, j_secondary, secondary_element_index)

    # Since the flux is computed in reference space, we do it separately for each element
    flux_primary = surface_flux(u_ll, u_rr_ll,
                                reference_normal_vector(primary_direction_index),
                                equations, cache.elements,
                                i_primary, j_primary, primary_element_index)
    flux_secondary = surface_flux(u_rr, u_ll_rr,
                                  reference_normal_vector(secondary_direction_index),
                                  equations, cache.elements,
                                  i_secondary, j_secondary, secondary_element_index)

    # Now we can update the surface flux values, where all vector variables are stored 
    # as contravariant components
    for v in eachvariable(equations)
        surface_flux_values[v, primary_node_index, primary_direction_index,
        primary_element_index] = flux_primary[v]
        surface_flux_values[v, secondary_node_index, secondary_direction_index,
        secondary_element_index] = flux_secondary[v]
    end
end

# This function passes the element container and node/element indices into the source term
function Trixi.calc_sources!(du, u, t, source_terms,
                             equations::AbstractCovariantEquations{2}, dg::DG, cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            u_local = Trixi.get_node_vars(u, equations, dg, i, j, element)
            x_local = Trixi.get_node_coords(cache.elements.node_coordinates, equations,
                                            dg,
                                            i, j, element)
            du_local = source_terms(u_local, x_local, t, equations,
                                    cache.elements, i, j, element)
            Trixi.add_to_node_vars!(du, du_local, equations, dg, i, j, element)
        end
    end

    return nothing
end

# Version for no sources
function Trixi.calc_sources!(du, u, t, source_terms::Nothing,
                             equations::AbstractCovariantEquations{2}, dg::DG, cache)
    return nothing
end

# Calculate time step based on maximum wave speeds for the covariant form
function Trixi.max_dt(u, t, mesh::P4estMesh{2}, constant_speed::False,
                      equations::AbstractCovariantEquations{2},
                      dg::DG, cache)

    # to avoid a division by zero if the speed vanishes everywhere,
    # e.g. for steady-state linear advection
    max_scaled_speed = nextfloat(zero(t))

    # Because the covariant form computes max_abs_speeds using the contravariant 
    # velocity components already, there is no need to transform them here
    for element in eachelement(dg, cache)
        max_lambda1 = max_lambda2 = zero(max_scaled_speed)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
            lambda1, lambda2 = Trixi.max_abs_speeds(u_node, equations, cache.elements,
                                                    i, j, element)
            max_lambda1 = max(max_lambda1, lambda1)
            max_lambda2 = max(max_lambda2, lambda2)
        end

        max_scaled_speed = max(max_scaled_speed, max_lambda1 + max_lambda2)
    end
    return 2 / (nnodes(dg) * max_scaled_speed)
end

function Trixi.calc_error_norms(func, u, t, analyzer, mesh::P4estMesh{2},
                                equations::AbstractCovariantEquations{2},
                                initial_condition, dg::DGSEM, cache, cache_analysis)
    (; weights) = dg.basis
    (; node_coordinates) = cache.elements

    # Set up data structures
    l2_error = zero(func(Trixi.get_node_vars(u, equations, dg, 1, 1, 1), equations))
    linf_error = copy(l2_error)
    total_volume = zero(real(mesh))

    # Iterate over all elements for error calculations
    for element in eachelement(dg, cache)

        # Calculate errors at each volume quadrature node
        for j in eachnode(dg), i in eachnode(dg)
            x = Trixi.get_node_coords(node_coordinates, equations, dg, i, j, element)

            u_exact = spherical2contravariant(initial_condition(x, t, equations),
                                              equations, cache.elements, i, j, element)

            u_numerical = Trixi.get_node_vars(u, equations, dg, i, j, element)

            diff = func(u_exact, equations) - func(u_numerical, equations)

            J = volume_element(cache.elements, i, j, element)

            l2_error += diff .^ 2 * (weights[i] * weights[j] * J)
            linf_error = @. max(linf_error, abs(diff))
            total_volume += weights[i] * weights[j] * J
        end
    end

    # For L2 error, divide by total volume
    l2_error = @. sqrt(l2_error / total_volume)

    return l2_error, linf_error
end
end # @muladd
