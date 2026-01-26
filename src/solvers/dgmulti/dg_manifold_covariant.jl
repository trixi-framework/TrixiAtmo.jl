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

# uses quadrature + projection to compute source terms.
function Trixi.calc_sources!(du, u, t, source_terms,
                             mesh, equations::AbstractCovariantEquations, dg::DGMulti, cache)
    rd = dg.basis
    md = mesh.md
    @unpack Pq = rd
    @unpack u_values, aux_quad_values, local_values_threaded = cache
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
                             mesh, equations::AbstractCovariantEquations, dg::DGMulti, cache)
    nothing
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

    # interpolate to quadrature points
    Trixi.apply_to_each_field(Trixi.mul_by!(rd.Vq), u_values, u)

    Trixi.@threaded for e in Trixi.eachelement(mesh, dg, cache)
        flux_values = local_values_threaded[Threads.threadid()]
        for i in 1:NDIMS
            for j in Trixi.eachindex(flux_values)
                u_node = u_values[j, e]
                aux_node = aux_quad_values[j, e]
                area_elem = area_element(aux_node, equations)
                flux_values[j] = flux(u_node, aux_node, i, equations)
            end

            Trixi.apply_to_each_field(Trixi.mul_by_accum!(weak_differentiation_matrices[i]),
                                      view(du, :, e), flux_values)
        end
    end
end

function Trixi.calc_interface_flux!(cache, surface_integral::SurfaceIntegralWeakForm,
                                    mesh::DGMultiMesh,
                                    have_nonconservative_terms::False,
                                    equations::AbstractCovariantEquations{NDIMS},
                                    dg::DGMulti{NDIMS_AMBIENT, <:Tri}) where {NDIMS_AMBIENT, NDIMS}
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
        # in reference element coordinates.
        ref_index = mod(face_node_index - 1, Nfq) + 1
        normal = SVector(ntuple(k -> nrstJ[k][ref_index], NDIMS))

        flux_face_values[idM] = surface_flux(uM, uP_transformed_to_M, auxM, auxM, normal, equations)
    end
end
