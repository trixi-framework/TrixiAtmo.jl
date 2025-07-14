# Compute coefficients for an initial condition that uses auxiliary variables
function Trixi.compute_coefficients!(::Nothing, u, initial_condition, t,
                                     mesh::DGMultiMesh, equations::AbstractCovariantEquations,
                                     dg::DGMulti, cache)
    md = mesh.md
    rd = dg.basis
    @unpack u_values, aux_quad_values = cache
    # evaluate the initial condition at quadrature points
    Trixi.@threaded for i in Trixi.each_quad_node_global(mesh, dg, cache)
        x_node = SVector(getindex.(md.xyzq, i))
        aux_node = aux_quad_values[i]
        u_values[i] = initial_condition(x_node, t, aux_node, equations)
    end

    # multiplying by Pq computes the L2 projection
    Trixi.apply_to_each_field(Trixi.mul_by!(rd.Pq), u, u_values)
end

# version for covariant equations on DGMultiMeshes
function Trixi.calc_volume_integral!(du, u, mesh::DGMultiMesh{NDIMS_AMBIENT, <:Trixi.NonAffine},
                               have_nonconservative_terms::False,
                               equations::AbstractCovariantEquations{NDIMS},
                               volume_integral::VolumeIntegralWeakForm, dg::DGMulti,
                               cache) where {NDIMS_AMBIENT, NDIMS}
    rd = dg.basis
    md = mesh.md
    (; weak_differentiation_matrices, u_values, aux_quad_values, local_values_threaded) = cache

    Jq = rd.Vq * md.J

    # interpolate to quadrature points
    Trixi.apply_to_each_field(Trixi.mul_by!(rd.Vq), u_values, u)

    Trixi.@threaded for e in Trixi.eachelement(mesh, dg, cache)
        flux_values = local_values_threaded[Threads.threadid()]
        for i in 1:NDIMS
            for j in Trixi.eachindex(flux_values)
                u_node = u_values[j, e]
                aux_node = aux_quad_values[j, e]
                area_elem = area_element(aux_node, equations)
                J_node = Jq[j, e]
                # Rescale the flux, such that the volume integral becomes independent of the thickness
                # of the spherical shell. We compute that thickness by taking the ratio of the element's
                # volume and the area element. 
                flux_values[j] = flux(u_node, aux_node, i, equations) * (J_node / area_elem)
            end

            Trixi.apply_to_each_field(Trixi.mul_by_accum!(weak_differentiation_matrices[i]),
                                      view(du, :, e), flux_values)
        end
    end
end

# Calculate flux at interface by passing auxiliary variables to the surface flux function
function Trixi.calc_interface_flux!(cache, surface_integral::SurfaceIntegralWeakForm,
                                    mesh::DGMultiMesh,
                                    have_nonconservative_terms::False,
                                    equations::AbstractCovariantEquations{NDIMS},
                                    dg::DGMulti{NDIMS_AMBIENT, <:Wedge}) where {NDIMS_AMBIENT, NDIMS}
    @unpack surface_flux = surface_integral
    md = mesh.md
    rd = dg.basis
    @unpack mapM, mapP, Jf = md
    @unpack nrstJ, Nfq = rd
    @unpack u_face_values, flux_face_values, aux_face_values = cache
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
        # in reference element coordinates. We only take the first NDIMS components, projecting
        # the normal vector onto the (reference) tangent space of the manifold, possibly yielding a zero vector.
        ref_index = mod(face_node_index - 1, Nfq) + 1
        normal = SVector(ntuple(k -> nrstJ[k][ref_index], NDIMS))
        normal = normal / norm(normal)
        
        # If the unprojected normal vector belonged one of the triangular faces of the reference wedge,
        # the normalized projected normal vector consists of NaNs. In this case, we skip the flux computation,
        if any(isnan, normal)
            continue
        end

        # Compute the norm of the normal vector transformed to physical coordinates
        basis = basis_contravariant(auxM, equations)
        norm_normal_transformed = norm(transpose(basis) * normal)
        area_elem = area_element(auxM, equations)

        # Scale the fluxes by inverse norm of transformed normal vector and
        # the area element to get the correct fluxes per unit area,
        # making the surface integral independent of the thickness of the spherical shell.
        factor = 1 / norm_normal_transformed * Jf[idM] / area_elem
        flux_face_values[idM] = surface_flux(uM, uP_transformed_to_M, auxM, auxM, normal, equations) * factor
    end
end
