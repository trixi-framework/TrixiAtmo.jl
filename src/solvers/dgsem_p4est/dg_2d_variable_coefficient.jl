@muladd begin

function Trixi.rhs!(du, u, t,
              mesh::Union{TreeMesh{2}, P4estMesh{2}},
              equations::AbstractVariableCoefficientEquations,
              boundary_conditions, source_terms::Source,
              dg::DG, cache) where {Source}
    # Reset du
    Trixi.@trixi_timeit Trixi.timer() "reset ∂u/∂t" Trixi.set_zero!(du, dg, cache)

    # Calculate volume integral
    Trixi.@trixi_timeit Trixi.timer() "volume integral" begin
        Trixi.calc_volume_integral!(du, u, mesh, equations,
                              dg.volume_integral, dg, cache)
    end

    # Prolong solution to interfaces
    Trixi.@trixi_timeit Trixi.timer() "prolong2interfaces" begin
        Trixi.prolong2interfaces!(cache, u, mesh, equations, dg)
    end

    # Calculate interface fluxes
    Trixi.@trixi_timeit Trixi.timer() "interface flux" begin
        Trixi.calc_interface_flux!(cache.elements.surface_flux_values, mesh,
                             Trixi.have_nonconservative_terms(equations),
                             have_aux_node_vars(equations), equations,
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

    # Prolong solution to mortars
    Trixi.@trixi_timeit Trixi.timer() "prolong2mortars" begin
        Trixi.prolong2mortars!(cache, u, mesh, equations,
                         dg.mortar, dg)
    end

    # Calculate mortar fluxes
    Trixi.@trixi_timeit Trixi.timer() "mortar flux" begin
        Trixi.calc_mortar_flux!(cache.elements.surface_flux_values, mesh,
                          Trixi.have_nonconservative_terms(equations),
                          have_aux_node_vars(equations), equations,
                          dg.mortar, dg.surface_integral, dg, cache)
    end

    # Calculate surface integrals
    Trixi.@trixi_timeit Trixi.timer() "surface integral" begin
        Trixi.calc_surface_integral!(du, u, mesh, equations,
                               dg.surface_integral, dg, cache)
    end

    # Apply Jacobian from mapping to reference element
    Trixi.@trixi_timeit Trixi.timer() "Jacobian" Trixi.apply_jacobian!(du, mesh, equations, dg, cache)

    # Calculate source terms
    Trixi.@trixi_timeit Trixi.timer() "source terms" begin
        Trixi.calc_sources!(du, u, t, source_terms, have_aux_node_vars(equations),
                      equations, dg, cache)
    end

    return nothing
end

function Trixi.calc_volume_integral!(du, u, mesh, equations,
                               volume_integral::VolumeIntegralWeakForm,
                               dg::DGSEM, cache)
    Trixi.@threaded for element in Trixi.eachelement(dg, cache)
        Trixi.weak_form_kernel!(du, u, element, mesh,
                          Trixi.have_nonconservative_terms(equations),
                          have_aux_node_vars(equations), equations,
                          dg, cache)
    end

    return nothing
end

@inline function Trixi.weak_form_kernel!(du, u,
                                         element,
                                         mesh::P4estMesh{2},
                                         nonconservative_terms::Trixi.False,
                                         have_aux_node_vars::True,
                                         equations::AbstractVariableCoefficientEquations{2},
                                         dg::DGSEM, cache, alpha = true)
    (; derivative_hat) = dg.basis
    (; contravariant_vectors) = cache.elements
    (; aux_node_vars) = cache.aux_vars
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
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_hat[ii, i],
                                            contravariant_flux1, equations, dg, ii, j,
                                            element)
        end

        # Compute the contravariant flux by taking the scalar product of the
        # second contravariant vector Ja^2 and the flux vector
        Ja21, Ja22 = Trixi.get_contravariant_vector(2, contravariant_vectors, i, j, element)
        contravariant_flux2 = Ja21 * flux1 + Ja22 * flux2
        for jj in eachnode(dg)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_hat[jj, j],
                                            contravariant_flux2, equations, dg, i, jj,
                                            element)
        end
    end
end 


function Trixi.calc_interface_flux!(surface_flux_values,
                              mesh::P4estMesh{2},
                              nonconservative_terms, have_aux_node_vars,
                              equations, surface_integral, dg::DG, cache)
    @unpack neighbor_ids, node_indices = cache.interfaces
    @unpack contravariant_vectors = cache.elements
    index_range = eachnode(dg)
    index_end = last(index_range)

    Trixi.@threaded for interface in Trixi.eachinterface(dg, cache)
        # Get element and side index information on the primary element
        primary_element = neighbor_ids[1, interface]
        primary_indices = node_indices[1, interface]
        primary_direction = Trixi.indices2direction(primary_indices)

        # Create the local i,j indexing on the primary element used to pull normal direction information
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
            # Get the normal direction on the primary element.
            # Contravariant vectors at interfaces in negative coordinate direction
            # are pointing inwards. This is handled by `get_normal_direction`.
            normal_direction = Trixi.get_normal_direction(primary_direction,
                                                    contravariant_vectors,
                                                    i_primary, j_primary,
                                                    primary_element)

            Trixi.calc_interface_flux!(surface_flux_values, mesh, nonconservative_terms,
                                 equations,
                                 surface_integral, dg, cache,
                                 interface, normal_direction,
                                 node, primary_direction, primary_element,
                                 node_secondary, secondary_direction, secondary_element)

            # Increment primary element indices to pull the normal direction
            i_primary += i_primary_step
            j_primary += j_primary_step
            # Increment the surface node index along the secondary element
            node_secondary += node_secondary_step
        end
    end

    return nothing
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
    (; aux_surface_node_vars) = cache.aux_vars
    (; surface_flux) = surface_integral

    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, primary_node_index,
                                            interface_index)

    # Get auxiliary variables
    aux_vars_ll, aux_vars_rr = get_surface_node_aux_vars(aux_surface_node_vars,
                                                        equations,
                                                        dg, primary_node_index,
                                                        interface_index)

    #compute flux 
    flux_ = surface_flux(u_ll, u_rr, aux_vars_ll, aux_vars_rr,
                        normal_direction, equations)

    # Store flux values
    for v in eachvariable(equations)
        surface_flux_values[v, primary_node_index, primary_direction_index, primary_element_index] = flux_[v]
        surface_flux_values[v, secondary_node_index, secondary_direction_index, secondary_element_index] = -flux_[v]
    end
end

function Trixi.prolong2boundaries!(cache, u,
                             mesh::P4estMesh{2},
                             equations, surface_integral, dg::DG)
    @unpack boundaries = cache
    index_range = eachnode(dg)

    Trixi.@threaded for boundary in Trixi.eachboundary(dg, cache)
        # Copy solution data from the element using "delayed indexing" with
        # a start value and a step size to get the correct face and orientation.
        element = boundaries.neighbor_ids[boundary]
        node_indices = boundaries.node_indices[boundary]

        i_node_start, i_node_step = Trixi.index_to_start_step_2d(node_indices[1], index_range)
        j_node_start, j_node_step = Trixi.index_to_start_step_2d(node_indices[2], index_range)

        i_node = i_node_start
        j_node = j_node_start
        for i in eachnode(dg)
            for v in eachvariable(equations)
                boundaries.u[v, i, boundary] = u[v, i_node, j_node, element]
            end
            i_node += i_node_step
            j_node += j_node_step
        end
    end

    return nothing
end


@inline function Trixi.calc_boundary_flux!(surface_flux_values, t, boundary_condition,#::Union{NamedTuple{(:x_neg, :x_pos, :y_neg, :y_pos)}, Dict{Symbol, Any}},
                                    mesh::P4estMesh{2},
                                    nonconservative_terms::Trixi.False,
                                    equations::AbstractVariableCoefficientEquations{2},
                                    surface_integral, dg::DG, cache,
                                    i_index, j_index,
                                    node_index, direction_index, element_index,
                                    boundary_index)
    (; boundaries) = cache
    (; node_coordinates, contravariant_vectors) = cache.elements
    (; surface_flux) = surface_integral
    (; aux_node_vars) = cache.aux_vars
    #@show stacktrace()
    # Extract solution data from boundary container
    u_inner = Trixi.get_node_vars(boundaries.u, equations, dg, node_index, boundary_index)
    aux_vars_inner = get_node_aux_vars(aux_node_vars, equations, dg, i_index, j_index, element_index)

    # Outward-pointing normal direction (not normalized)
    normal_direction = Trixi.get_normal_direction(direction_index, contravariant_vectors,
                                            i_index, j_index, element_index)

    # Coordinates at boundary node
    x = Trixi.get_node_coords(node_coordinates, equations, dg, i_index, j_index,
                        element_index)
    flux_ = boundary_condition(u_inner, aux_vars_inner, normal_direction, x, t,
                            surface_flux, equations)

    # Copy flux to element storage in the correct orientation
    for v in eachvariable(equations)
        surface_flux_values[v, node_index, direction_index, element_index] = flux_[v]
    end
end

#function Trixi.calc_boundary_flux!(cache, t,
#                                   boundary_conditions::NamedTuple{(:x_neg, :x_pos, :y_neg, :y_pos)},
#                                   mesh::P4estMesh{2},
#                                   equations::AbstractVariableCoefficientEquations{2},
#                                   surface_integral, dg::DG, cache_extra...)
#
#    # Call your auxiliary-variable version here
#    #Trixi.calc_boundary_flux_with_aux!(cache, t, boundary_conditions, mesh, equations, surface_integral, dg, cache_extra...)
#end

function Trixi.calc_mortar_flux!(surface_flux_values,
                           mesh::P4estMesh{2},
                           nonconservative_terms, 
                           have_aux_node_vars::True, equations,
                           mortar_l2::Trixi.LobattoLegendreMortarL2,
                           surface_integral, dg::DG, cache)
    @unpack neighbor_ids, node_indices = cache.mortars
    @unpack contravariant_vectors = cache.elements
    @unpack (fstar_primary_upper_threaded, fstar_primary_lower_threaded,
    fstar_secondary_upper_threaded, fstar_secondary_lower_threaded) = cache
    @unpack aux_mortar_node_vars = cache.aux_vars

    index_range = eachnode(dg)

    Trixi.@threaded for mortar in Trixi.eachmortar(dg, cache)
        # Choose thread-specific pre-allocated container
        fstar_primary = (fstar_primary_lower_threaded[Threads.threadid()],
                         fstar_primary_upper_threaded[Threads.threadid()])

        fstar_secondary = (fstar_secondary_lower_threaded[Threads.threadid()],
                           fstar_secondary_upper_threaded[Threads.threadid()])

        # Get index information on the small elements
        small_indices = node_indices[1, mortar]
        small_direction = Trixi.indices2direction(small_indices)

        i_small_start, i_small_step = Trixi.index_to_start_step_2d(small_indices[1],
                                                             index_range)
        j_small_start, j_small_step = Trixi.index_to_start_step_2d(small_indices[2],
                                                             index_range)

        for position in 1:2
            i_small = i_small_start
            j_small = j_small_start
            element = neighbor_ids[position, mortar]
            for node in eachnode(dg)
                # Get the normal direction on the small element.
                # Note, contravariant vectors at interfaces in negative coordinate direction
                # are pointing inwards. This is handled by `get_normal_direction`.
                normal_direction = Trixi.get_normal_direction(small_direction,
                                                        contravariant_vectors,
                                                        i_small, j_small, element)

                Trixi.calc_mortar_flux!(fstar_primary, fstar_secondary,
                                  mesh, nonconservative_terms, equations,
                                  surface_integral, dg, cache,
                                  mortar, position, normal_direction,
                                  node)

                i_small += i_small_step
                j_small += j_small_step
            end
        end

        # Buffer to interpolate flux values of the large element to before
        # copying in the correct orientation
        u_buffer = cache.u_threaded[Threads.threadid()]

        # in calc_interface_flux!, the interface flux is computed once over each
        # interface using the normal from the "primary" element. The result is then
        # passed back to the "secondary" element, flipping the sign to account for the
        # change in the normal direction. For mortars, this sign flip occurs in
        # "mortar_fluxes_to_elements!" instead.
        Trixi.mortar_fluxes_to_elements!(surface_flux_values,
                                   mesh, equations, mortar_l2, dg, cache,
                                   mortar, fstar_primary, fstar_secondary,
                                   u_buffer)
    end

    return nothing
end

# Inlined version of the mortar flux computation on small elements for conservation laws
@inline function Trixi.calc_mortar_flux!(fstar_primary, fstar_secondary,
                                   mesh::P4estMesh{2},
                                   nonconservative_terms::False, equations,
                                   surface_integral, dg::DG, cache,
                                   mortar_index, position_index, normal_direction,
                                   node_index)
    @unpack u = cache.mortars
    @unpack aux = cache.aux_vars
    @unpack surface_flux = surface_integral

    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, position_index,
                                       node_index, mortar_index)
    aux_ll, aux_rr = get_aux_surface_node_vars(aux, equations, dg,
                                                    position_index,
                                       node_index, mortar_index)

    flux = surface_flux(u_ll, u_rr, aux_ll, aux_rr, normal_direction, equations)

    # Copy flux to buffer
    Trixi.set_node_vars!(fstar_primary[position_index], flux, equations, dg, node_index)
    Trixi.set_node_vars!(fstar_secondary[position_index], flux, equations, dg, node_index)
end

function Trixi.calc_sources!(du, u, t, source_terms::Nothing, have_aux_node_vars::True,
                       equations::AbstractEquations{2}, dg::DG, cache)
    return nothing
end
function Trixi.calc_sources!(du, u, t, source_terms, have_aux_node_vars::True,
                       equations::AbstractEquations{2}, dg::DG, cache)
    @unpack node_coordinates = cache.elements
    @unpack aux_node_vars = cache.aux_vars

    Trixi.@threaded for element in Trixi.eachelement(dg, cache)
        for j in Trixi.eachnode(dg), i in Trixi.eachnode(dg)
            u_local = Trixi.get_node_vars(u, equations, dg, i, j, element)
            aux_local = get_aux_node_vars(aux_node_vars, equations, dg,
                                          i, j, element)
            x_local = Trixi.get_node_coords(node_coordinates, equations, dg,
                                      i, j, element)
            du_local = source_terms(u_local, aux_local, x_local, t, equations)
            Trixi.add_to_node_vars!(du, du_local, equations, dg, i, j, element)
        end
    end

    return nothing
end
end 