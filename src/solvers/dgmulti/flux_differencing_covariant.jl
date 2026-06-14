@muladd begin
#! format: noindent

function Trixi.create_cache(mesh::DGMultiMesh, equations::AbstractCovariantEquations,
                            dg::Trixi.DGMultiFluxDiffSBP,
                            RealT, metric_terms, auxiliary_field,
                            uEltype)
    rd = dg.basis
    md = mesh.md

    # for use with flux differencing schemes
    Qrst_skew = Trixi.compute_flux_differencing_SBP_matrices(dg)

    lift_scalings = rd.wf ./ rd.wq[rd.Fmask] # lift scalings for diag-norm SBP operators

    nvars = nvariables(equations)
    naux = n_aux_node_vars(equations)
    # Use an array of SVectors (chunks of `nvars` are contiguous in memory) to speed up flux differencing
    du_local_threaded = [zeros(SVector{nvars, uEltype}, rd.Nq)
                         for _ in 1:Threads.maxthreadid()]

    solution_container = Trixi.initialize_dgmulti_solution_container(mesh, equations,
                                                                     dg,
                                                                     uEltype)

    aux_values = Trixi.allocate_nested_array(uEltype, naux, size(md.x), dg)
    aux_quad_values = Trixi.allocate_nested_array(uEltype, naux, size(md.xq), dg)
    aux_face_values = Trixi.allocate_nested_array(uEltype, naux, size(md.xf), dg)
    # use StructArray for better memory access in flux differencing
    aux_values, aux_quad_values, aux_face_values = StructArray.((aux_values,
                                                                 aux_quad_values,
                                                                 aux_face_values))
    auxiliary_container = (; aux_values, aux_quad_values, aux_face_values)

    init_auxiliary_node_variables!(aux_values, mesh, equations, dg, metric_terms,
                                   auxiliary_field)
    aux_quad_values .= aux_values
    Trixi.apply_to_each_field(Trixi.mul_by!(rd.Vf), aux_face_values, aux_values)

    invJ = inv.(area_element.(aux_quad_values, equations))

    return (; md, Qrst_skew,
            invJ = invJ, lift_scalings, inv_wq = inv.(rd.wq),
            solution_container, auxiliary_container,
            du_local_threaded,)
end

@inline function Trixi.local_flux_differencing!(du_local, u_local, aux_local,
                                                element_index,
                                                have_nonconservative_terms::True,
                                                volume_flux,
                                                has_sparse_operators::False, mesh,
                                                equations::AbstractCovariantEquations{NDIMS},
                                                dg, cache) where {NDIMS}
    @unpack Qrst_skew = cache
    flux_conservative, flux_nonconservative = volume_flux
    row_ids = axes(first(Qrst_skew), 1)
    for i in row_ids
        u_i = u_local[i]
        aux_i = aux_local[i]
        for j in row_ids
            normal_direction = SVector(ntuple(d -> Qrst_skew[d][i, j], Val(NDIMS)))
            # We use the symmetry of the volume flux and the anti-symmetry
            # of the derivative operator to save half of the volume flux
            # computations.
            if j > i
                u_j = u_local[j]
                aux_j = aux_local[j]
                AF_ij = 2 * flux_conservative(u_i, u_j, aux_i, aux_j, normal_direction,
                                          equations)
                du_local[i] = du_local[i] + AF_ij
                du_local[j] = du_local[j] - AF_ij # Due to skew-symmetry
            end
            # Non-conservative terms use the full (non-symmetric) loop.
            # The 0.5f0 factor on the normal direction is necessary for the nonconservative 
            # fluxes based on the interpretation of global SBP operators.  
            # See also `calc_interface_flux!` with `have_nonconservative_terms::True` 
            # in src/solvers/dgsem_tree/dg_1d.jl
            f_nc = flux_nonconservative(u_i, u_local[j], aux_i, aux_local[j],
                                        0.5f0 * normal_direction,
                                        equations)
            du_local[i] = du_local[i] + 2 * f_nc
        end
    end
end

@inline function Trixi.volume_integral_kernel!(du, u, element, mesh::DGMultiMesh,
                                               have_nonconservative_terms,
                                               equations::AbstractCovariantEquations,
                                               volume_integral::VolumeIntegralFluxDifferencing,
                                               dg::Trixi.DGMultiFluxDiffSBP, cache,
                                               alpha = true)
    @unpack du_local_threaded, inv_wq = cache
    (; aux_values) = cache.auxiliary_container

    du_local = du_local_threaded[Threads.threadid()]
    fill!(du_local, zero(eltype(du_local)))
    u_local = view(u, :, element)
    aux_local = view(aux_values, :, element)

    Trixi.local_flux_differencing!(du_local, u_local, aux_local, element,
                                   have_nonconservative_terms,
                                   volume_integral.volume_flux,
                                   Trixi.has_sparse_operators(dg),
                                   mesh, equations, dg, cache)

    for i in Trixi.each_quad_node(mesh, dg, cache)
        du[i, element] = du[i, element] + alpha * du_local[i] * inv_wq[i]
    end

    return nothing
end

# Specialize since `u_values` isn't computed for DGMultiFluxDiffSBP solvers.
function Trixi.calc_sources!(du, u, t, source_terms,
                             mesh, equations::AbstractCovariantEquations,
                             dg::Trixi.DGMultiFluxDiffSBP, cache)
    md = mesh.md
    (; aux_values) = cache.auxiliary_container

    Trixi.@threaded for e in Trixi.eachelement(mesh, dg, cache)
        for i in Trixi.each_quad_node(mesh, dg, cache)
            du[i, e] += source_terms(u[i, e], SVector(getindex.(md.xyzq, i, e)), t,
                                     aux_values[i, e], equations)
        end
    end
end
end # @muladd
