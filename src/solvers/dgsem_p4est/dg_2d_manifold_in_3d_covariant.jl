@muladd begin
#! format: noindent

"""
    Get the normal vector to the reference quadrilateral at a given node
"""
@inline function reference_normal_vector(direction)
    orientation = (direction + 1) >> 1
    sign = isodd(direction) ? -1 : 1

    if orientation == 1
        return SVector(sign, 0.0f0, 0.0f0)
    else
        return SVector(0.0f0, sign, 0.0f0)
    end
end

"""
   Get Cartesian coordinates at a point on the manifold
"""
@inline function Trixi.get_node_coords(x, ::AbstractCovariantEquations2D, ::DG,
                                       indices...)
    return SVector(ntuple(@inline(idx->x[idx, indices...]), 3))
end

"""
    Evaluate the state with Cartesian vector components, then transform to contravariant components for use within the solver.
"""
function Trixi.compute_coefficients!(u, func, t, mesh::Trixi.AbstractMesh{2},
                                     equations::AbstractCovariantEquations2D, dg::DG,
                                     cache)
    for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            x_node = view(cache.elements.node_coordinates, :, i, j, element)
            u_car = func(x_node, t, equations)
            u_node = cartesian2contravariant(u_car, equations, i, j, element, cache)
            Trixi.set_node_vars!(u, u_node, equations, dg, i, j, element)
        end
    end
end

"""
    Weak form kernel for fluxes already in contravariant components
"""
@inline function Trixi.weak_form_kernel!(du, u, element,
                                         mesh::Union{P4estMesh{2},
                                                     StructuredMesh{2}, T8codeMesh{2},
                                                     UnstructuredMesh2D},
                                         nonconservative_terms::False,
                                         equations::AbstractCovariantEquations2D,
                                         dg::DGSEM, cache, alpha = true)
    (; derivative_dhat) = dg.basis
    (; inverse_jacobian) = cache.elements

    for j in eachnode(dg), i in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        J_node = 1 / inverse_jacobian[i, j, element]

        contravariant_flux1 = J_node * flux(u_node, 1, equations)
        contravariant_flux2 = J_node * flux(u_node, 2, equations)

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

"""
    Weak form kernel for fluxes already in contravariant components
"""
@inline function Trixi.flux_differencing_kernel!(du, u,
                                                 element,
                                                 mesh::Union{StructuredMesh{2},
                                                            StructuredMeshView{2},
                                                            UnstructuredMesh2D, P4estMesh{2},
                                                            T8codeMesh{2}},
                                                 nonconservative_terms::False, 
                                                 equations::AbstractCovariantEquations2D,
                                                 volume_flux, dg::DGSEM, cache, 
                                                 alpha = true)
    (; derivative_split) = dg.basis
    (; inverse_jacobian) = cache.elements

    for j in eachnode(dg), i in eachnode(dg)
        u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
        J_node = 1 / inverse_jacobian[i, j, element]

        # x direction
        for ii in (i + 1):nnodes(dg)
            u_node_ii = Trixi.get_node_vars(u, equations, dg, ii, j, element)
            J_avg = 0.5 * (J_node +  1 / inverse_jacobian[ii, j, element])
            flux1 = J_avg * volume_flux(u_node, u_node_ii, 1, equations)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[i, ii], flux1,
                                       equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[ii, i], flux1,
                                       equations, dg, ii, j, element)
        end

        # y direction
        for jj in (j + 1):nnodes(dg)
            u_node_jj = Trixi.get_node_vars(u, equations, dg, i, jj, element)
            J_avg = 0.5 * (J_node +  1 / inverse_jacobian[i, jj, element])
            flux2 = J_avg * volume_flux(u_node, u_node_jj, 2, equations)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[j, jj], flux2,
                                       equations, dg, i, j, element)
            Trixi.multiply_add_to_node_vars!(du, alpha * derivative_split[jj, j], flux2,
                                       equations, dg, i, jj, element)
        end
    end

    return nothing
end

"""
    Interface flux for the covariant form, where the numerical flux is computed in the direction of the 2D reference normal
"""
function Trixi.calc_interface_flux!(surface_flux_values,
                                    mesh::P4estMesh{2}, nonconservative_terms,
                                    equations::AbstractCovariantEquations2D,
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

"""
    Pointwise computation of the interface flux for the covariant form
"""
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
    (; inverse_jacobian) = cache.elements
    (; u) = cache.interfaces
    (; surface_flux) = surface_integral

    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, primary_node_index,
                                             interface_index)

    # Convert to contravariant components on each element
    u_ll_car = contravariant2cartesian(u_ll, equations, i_primary, j_primary,
                                       primary_element_index, cache)
    u_rr_car = contravariant2cartesian(u_rr, equations, i_secondary, j_secondary,
                                       secondary_element_index, cache)
    u_rr_ll = cartesian2contravariant(u_rr_car, equations, i_primary, j_primary,
                                      primary_element_index, cache)
    u_ll_rr = cartesian2contravariant(u_ll_car, equations, i_secondary, j_secondary,
                                      secondary_element_index, cache)

    # Get the normal vector in reference space to each element
    reference_normal_primary = reference_normal_vector(primary_direction_index)
    reference_normal_secondary = reference_normal_vector(secondary_direction_index)

    # Since the flux is computed in reference space, we do it separately for each element
    flux_primary = surface_flux(u_ll, u_rr_ll, reference_normal_primary, equations) /
                   inverse_jacobian[i_primary, j_primary, primary_element_index]
    flux_secondary = surface_flux(u_rr, u_ll_rr, reference_normal_secondary,
                                  equations) /
                     inverse_jacobian[i_secondary, j_secondary, secondary_element_index]

    for v in eachvariable(equations)
        surface_flux_values[v, primary_node_index, primary_direction_index,
        primary_element_index] = flux_primary[v]
        surface_flux_values[v, secondary_node_index, secondary_direction_index,
        secondary_element_index] = flux_secondary[v]
    end
end

"""
    Max time step for wave speeds given directly in contravariant components
"""
function Trixi.max_dt(u, t,
                      mesh::Union{P4estMesh{2}, StructuredMesh{2}, T8codeMesh{2},
                                  UnstructuredMesh2D},
                      constant_speed::False,
                      equations::AbstractCovariantEquations2D, dg::DG, cache)

    # to avoid a division by zero if the speed vanishes everywhere,
    # e.g. for steady-state linear advection
    max_scaled_speed = nextfloat(zero(t))

    # Because the covariant form computes max_abs_speeds using the contravariant 
    # velocity components already, there is no need to transform them here
    for element in eachelement(dg, cache)
        max_lambda1 = max_lambda2 = zero(max_scaled_speed)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = Trixi.get_node_vars(u, equations, dg, i, j, element)
            lambda1, lambda2 = Trixi.max_abs_speeds(u_node, equations)

            max_lambda1 = max(max_lambda1, lambda1)
            max_lambda2 = max(max_lambda2, lambda2)
        end

        max_scaled_speed = max(max_scaled_speed, max_lambda1 + max_lambda2)
    end

    return 2 / (nnodes(dg) * max_scaled_speed)
end

"""
    Specialize create_cache_analysis for 3D vector of node coordinates with a 2D equation
"""
function Trixi.create_cache_analysis(analyzer,
                                     mesh::Union{StructuredMesh{2}, UnstructuredMesh2D,
                                                 P4estMesh{2}, T8codeMesh{2}},
                                     equations::AbstractCovariantEquations2D, dg::DG,
                                     cache,
                                     RealT, uEltype)

    # pre-allocate buffers
    # We use `StrideArray`s here since these buffers are used in performance-critical
    # places and the additional information passed to the compiler makes them faster
    # than native `Array`s.
    u_local = StrideArray(undef, uEltype,
                          StaticInt(nvariables(equations)), StaticInt(nnodes(analyzer)),
                          StaticInt(nnodes(analyzer)))
    u_tmp1 = StrideArray(undef, uEltype,
                         StaticInt(nvariables(equations)), StaticInt(nnodes(analyzer)),
                         StaticInt(nnodes(dg)))
    x_local = StrideArray(undef, RealT,
                          StaticInt(3), StaticInt(nnodes(analyzer)),
                          StaticInt(nnodes(analyzer)))
    x_tmp1 = StrideArray(undef, RealT,
                         StaticInt(3), StaticInt(nnodes(analyzer)),
                         StaticInt(nnodes(dg)))
    jacobian_local = StrideArray(undef, RealT,
                                 StaticInt(nnodes(analyzer)),
                                 StaticInt(nnodes(analyzer)))
    jacobian_tmp1 = StrideArray(undef, RealT,
                                StaticInt(nnodes(analyzer)), StaticInt(nnodes(dg)))

    return (; u_local, u_tmp1, x_local, x_tmp1, jacobian_local, jacobian_tmp1)
end

"""
    Calculate error norms for systems on manifolds
"""
function Trixi.calc_error_norms(func, u, t, analyzer,
                                mesh::Union{StructuredMesh{2}, UnstructuredMesh2D,
                                            P4estMesh{2},
                                            T8codeMesh{2}},
                                equations::AbstractCovariantEquations2D,
                                initial_condition, dg::DGSEM, cache, cache_analysis)
    (; weights) = dg.basis
    (; node_coordinates, inverse_jacobian) = cache.elements

    # Set up data structures
    l2_error = zero(func(Trixi.get_node_vars(u, equations, dg, 1, 1, 1), equations))
    linf_error = copy(l2_error)
    total_volume = zero(real(mesh))

    # Iterate over all elements for error calculations
    for element in eachelement(dg, cache)
        # Calculate errors at each volume quadrature node
        for j in eachnode(dg), i in eachnode(dg)
            x = Trixi.get_node_coords(node_coordinates, equations, dg, i, j, element)
            u_exact = initial_condition(x, t, equations)

            diff = cartesian2contravariant(func(u_exact, equations),
                                           equations, i, j, element, cache) -
                   func(Trixi.get_node_vars(u, equations, dg, i, j, element), equations)

            J = inv(abs(inverse_jacobian[i, j, element]))

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
