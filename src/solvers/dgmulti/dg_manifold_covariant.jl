# Compute coefficients for an initial condition that uses auxiliary variables
function Trixi.compute_coefficients!(::Nothing, u, initial_condition, t,
                                     mesh::DGMultiMesh,
                                     equations::AbstractCovariantEquations,
                                     dg::DGMulti, cache)
    md = mesh.md
    rd = dg.basis
    (; u_values) = cache.solution_container
    (; aux_quad_values) = cache.auxiliary_container
    # evaluate the initial condition at quadrature points
    Trixi.@threaded for i in Trixi.each_quad_node_global(mesh, dg, cache)
        x_node = SVector(getindex.(md.xyzq, i))
        aux_node = aux_quad_values[i]
        u_values[i] = initial_condition(x_node, t, aux_node, equations)
    end

    # multiplying by Pq computes the L2 projection
    Trixi.apply_to_each_field(Trixi.mul_by!(rd.Pq), u, u_values)
end

# Calculate the source contribution, passing auxiliary variables to the source term function.
function Trixi.calc_sources!(du, u, t, source_terms,
                             mesh, equations::AbstractCovariantEquations, dg::DGMulti,
                             cache)
    rd = dg.basis
    md = mesh.md
    @unpack Pq = rd
    (; u_values, local_values_threaded) = cache.solution_container
    (; aux_quad_values) = cache.auxiliary_container
    Trixi.@threaded for e in Trixi.eachelement(mesh, dg, cache)
        source_values = local_values_threaded[Threads.threadid()]

        u_e = view(u_values, :, e) # u_values should already be computed from volume integral
        aux_e = view(aux_quad_values, :, e)

        for i in Trixi.each_quad_node(mesh, dg, cache)
            source_values[i] = source_terms(u_e[i], SVector(getindex.(md.xyzq, i, e)),
                                            t, aux_e[i], equations)
        end
        Trixi.apply_to_each_field(Trixi.mul_by_accum!(Pq), view(du, :, e), source_values)
    end

    return nothing
end

function Trixi.calc_sources!(du, u, t, source_term::Nothing,
                             mesh, equations::AbstractCovariantEquations, dg::DGMulti,
                             cache)
    nothing
end

# Affine mesh version for covariant equations - dispatcher
function Trixi.volume_integral_kernel!(du, u, element,
                                       mesh::DGMultiMesh{NDIMS_AMBIENT, <:Trixi.Affine},
                                       have_nonconservative_terms::False,
                                       equations::AbstractCovariantEquations{NDIMS},
                                       volume_integral::VolumeIntegralWeakForm, dg::DGMulti,
                                       cache) where {NDIMS_AMBIENT, NDIMS}
    volume_integral_covariant_kernel!(du, u, element, mesh, have_nonconservative_terms,
                                      equations, volume_integral, dg, cache)
end

# Non-affine mesh version for covariant equations - dispatcher
function Trixi.volume_integral_kernel!(du, u, element,
                                       mesh::DGMultiMesh{NDIMS_AMBIENT, <:Trixi.NonAffine},
                                       have_nonconservative_terms::False,
                                       equations::AbstractCovariantEquations{NDIMS},
                                       volume_integral::VolumeIntegralWeakForm, dg::DGMulti,
                                       cache) where {NDIMS_AMBIENT, NDIMS}
    volume_integral_covariant_kernel!(du, u, element, mesh, have_nonconservative_terms,
                                      equations, volume_integral, dg, cache)
end

# Volume integral kernel for covariant equations.
function volume_integral_covariant_kernel!(du, u, element,
                                           mesh::DGMultiMesh,
                                           have_nonconservative_terms::False,
                                           equations::AbstractCovariantEquations{NDIMS},
                                           volume_integral::VolumeIntegralWeakForm,
                                           dg::DGMulti,
                                           cache) where {NDIMS}
    (; weak_differentiation_matrices) = cache
    (; u_values, local_values_threaded) = cache.solution_container
    (; aux_quad_values) = cache.auxiliary_container

    flux_values = local_values_threaded[Threads.threadid()]
    for i in 1:NDIMS
        for j in Trixi.eachindex(flux_values)
            u_node = u_values[j, element]
            aux_node = aux_quad_values[j, element]
            flux_values[j] = flux(u_node, aux_node, i, equations)
        end

        Trixi.apply_to_each_field(Trixi.mul_by_accum!(weak_differentiation_matrices[i]),
                                  view(du, :, element), flux_values)
    end
end

function Trixi.calc_interface_flux!(cache,
                                    surface_integral::SurfaceIntegralWeakForm,
                                    mesh::DGMultiMesh,
                                    have_nonconservative_terms::False,
                                    equations::AbstractCovariantEquations{NDIMS},
                                    dg::DGMulti{NDIMS_AMBIENT, <:Tri}) where {NDIMS_AMBIENT,
                                                                              NDIMS}
    @unpack surface_flux = surface_integral
    md = mesh.md
    rd = dg.basis
    @unpack mapM, mapP, Jf = md
    @unpack nrstJ, Nfq = rd
    (; u_face_values, flux_face_values) = cache.solution_container
    (; aux_face_values) = cache.auxiliary_container
    Trixi.@threaded for face_node_index in Trixi.each_face_node_global(mesh, dg, cache)

        # inner (idM -> minus) and outer (idP -> plus) indices
        idM, idP = mapM[face_node_index], mapP[face_node_index]
        uM = u_face_values[idM]
        uP = u_face_values[idP]
        auxM = aux_face_values[idM]
        auxP = aux_face_values[idP]

        # Transform uP to the same coordinate system as uM
        uP_global = contravariant2global(uP, auxP, equations)
        uP_transformed_to_M = global2contravariant(uP_global, auxM, equations)

        # Compute ref_index from face_node_index in order to access covariant normal vector
        # in reference element coordinates.
        ref_index = mod(face_node_index - 1, Nfq) + 1
        normal = SVector{NDIMS}(getindex.(nrstJ, ref_index))

        flux_face_values[idM] = surface_flux(uM, uP_transformed_to_M, auxM, auxM, normal,
                                             equations)
    end
end

function Trixi.calc_interface_flux!(cache, surface_integral::SurfaceIntegralWeakForm,
                                    mesh::DGMultiMesh,
                                    have_nonconservative_terms::True,
                                    equations::AbstractCovariantEquations{NDIMS},
                                    dg::DGMulti{NDIMS}) where {NDIMS}
    flux_conservative, flux_nonconservative = surface_integral.surface_flux
    md = mesh.md
    rd = dg.basis
    @unpack mapM, mapP = md
    @unpack nrstJ, Nfq = rd
    (; u_face_values, flux_face_values) = cache.solution_container
    (; aux_face_values) = cache.auxiliary_container

    Trixi.@threaded for face_node_index in Trixi.each_face_node_global(mesh, dg, cache)

        # inner (idM -> minus) and outer (idP -> plus) indices
        idM, idP = mapM[face_node_index], mapP[face_node_index]
        uM = u_face_values[idM]
        auxM = aux_face_values[idM]

        # compute flux if node is not a boundary node
        if idM != idP
            uP = u_face_values[idP]
            auxP = aux_face_values[idP]

            # Transform uP to the same coordinate system as uM
            uP_global = contravariant2global(uP, auxP, equations)
            uP_transformed_to_M = global2contravariant(uP_global, auxM, equations)

            # Compute ref_index from face_node_index in order to access covariant normal vector
            # in reference element coordinates.
            ref_index = mod(face_node_index - 1, Nfq) + 1
            normal = SVector{NDIMS}(getindex.(nrstJ, ref_index))

            conservative_part = flux_conservative(uM, uP_transformed_to_M, auxM, auxM,
                                                  normal, equations)

            # Two notes on the use of `flux_nonconservative`:
            # 1. In contrast to other mesh types, only one nonconservative part needs to be
            #    computed since we loop over the elements, not the unique interfaces.
            nonconservative_part = flux_nonconservative(uM, uP_transformed_to_M, auxM, auxM,
                                                        normal, equations)
            # The factor 0.5 is necessary for the nonconservative fluxes based on the
            # interpretation of global SBP operators.
            flux_face_values[idM] = (conservative_part + 0.5 * nonconservative_part)
        end
    end

    return nothing
end

function Trixi.calc_single_boundary_flux!(cache, t, boundary_condition, boundary_key, mesh,
                                          have_nonconservative_terms::False,
                                          equations::AbstractCovariantEquations{NDIMS},
                                          dg::DGMulti) where {NDIMS}
    rd = dg.basis
    md = mesh.md
    (; u_face_values, flux_face_values) = cache.solution_container
    (; aux_face_values) = cache.auxiliary_container
    @unpack xyzf = md
    @unpack surface_flux = dg.surface_integral

    # reshape face/normal arrays to have size = (num_points_on_face, num_faces_total).
    # mesh.boundary_faces indexes into the columns of these face-reshaped arrays.
    num_faces = StartUpDG.num_faces(rd.element_type)
    num_pts_per_face = rd.Nfq ÷ num_faces
    num_faces_total = num_faces * md.num_elements

    # This function was originally defined as
    # `reshape_by_face(u) = reshape(view(u, :), num_pts_per_face, num_faces_total)`.
    # This results in allocations due to https://github.com/JuliaLang/julia/issues/36313.
    # To avoid allocations, we use Tim Holy's suggestion:
    # https://github.com/JuliaLang/julia/issues/36313#issuecomment-782336300.
    reshape_by_face(u) = Base.ReshapedArray(u, (num_pts_per_face, num_faces_total), ())

    u_face_values = reshape_by_face(u_face_values)
    aux_face_values = reshape_by_face(aux_face_values)
    flux_face_values = reshape_by_face(flux_face_values)
    xyzf = reshape_by_face.(xyzf)
    nrstJ = map(nrstJ -> Base.ReshapedArray(nrstJ, (num_pts_per_face, num_faces), ()),
                rd.nrstJ)

    # loop through boundary faces, which correspond to columns of reshaped u_face_values, ...
    for f in mesh.boundary_faces[boundary_key]
        for i in Base.OneTo(num_pts_per_face)
            ref_face = mod(f - 1, num_faces) + 1
            face_normal = SVector{NDIMS}(getindex.(nrstJ, i, ref_face))
            face_coordinates = SVector{NDIMS}(getindex.(xyzf, i, f))
            flux_face_values[i, f] = boundary_condition(u_face_values[i, f],
                                                        aux_face_values[i, f],
                                                        face_normal, face_coordinates,
                                                        t,
                                                        surface_flux, equations)
        end
    end

    # Note: modifying the values of the reshaped array modifies the values of cache.solution_container.flux_face_values.
    # However, we don't have to re-reshape, since cache.solution_container.flux_face_values still retains its original shape.
end
