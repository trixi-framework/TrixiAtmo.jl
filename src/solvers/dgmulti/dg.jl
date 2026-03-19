# Constructs cache variables including auxiliary variables for covariant equations and DGMultiMeshes
function Trixi.create_cache(mesh::DGMultiMesh{NDIMS}, equations::AbstractCovariantEquations,
                            dg::Trixi.DGMultiWeakForm,
                            RealT, metric_terms, auxiliary_field,
                            uEltype) where {NDIMS}
    rd = dg.basis
    md = mesh.md

    # volume quadrature weights, volume interpolation matrix, mass matrix, differentiation matrices
    @unpack wq, Vq, M, Drst = rd

    # ∫f(u) * dv/dx_i = ∑_j (Vq*Drst[i])'*diagm(wq)*(rstxyzJ[i,j].*f(Vq*u))
    weak_differentiation_matrices = map(D -> -M \ ((Vq * D)' * Diagonal(wq)), Drst)

    nvars = nvariables(equations)
    naux = n_aux_node_vars(equations)

    # We are duplicating the contents of solution_container in the top-level cache, but 
    # note that no actual data is being copied here, just references to the same arrays.
    u_values = Trixi.allocate_nested_array(uEltype, nvars, size(md.xq), dg)
    u_face_values = Trixi.allocate_nested_array(uEltype, nvars, size(md.xf), dg)
    flux_face_values = Trixi.allocate_nested_array(uEltype, nvars, size(md.xf), dg)
    local_values_threaded = [Trixi.allocate_nested_array(uEltype, nvars, (rd.Nq,), dg)
                             for _ in 1:Threads.maxthreadid()]
    solution_container = (; u_values, u_face_values, flux_face_values,
                          local_values_threaded)

    # To parallel the solution container, we create an auxiliary container.
    aux_values = Trixi.allocate_nested_array(uEltype, naux, size(md.x), dg)
    aux_quad_values = Trixi.allocate_nested_array(uEltype, naux, size(md.xq), dg)
    aux_face_values = Trixi.allocate_nested_array(uEltype, naux, size(md.xf), dg)
    auxiliary_container = (; aux_values, aux_quad_values, aux_face_values)

    if typeof(rd.approximation_type) <:
       Union{SBP, Trixi.AbstractNonperiodicDerivativeOperator}
        lift_scalings = rd.wf ./ rd.wq[rd.Fmask] # lift scalings for diag-norm SBP operators
    else
        lift_scalings = nothing
    end

    # For curved meshes, we interpolate geometric terms from nodal points to quadrature points.
    # For affine meshes, we just access one element of this interpolated data.
    dxidxhatj = map(x -> rd.Vq * x, md.rstxyzJ)

    init_auxiliary_node_variables!(aux_values, mesh, equations, dg, auxiliary_field)

    # Interpolate auxiliary variables to quadrature and face points
    Trixi.apply_to_each_field(Trixi.mul_by!(rd.Vq), aux_quad_values, aux_values)
    Trixi.apply_to_each_field(Trixi.mul_by!(rd.Vf), aux_face_values, aux_values)

    # interpolate J to quadrature points for weight-adjusted DG (WADG)
    invJ = inv.(area_element.(aux_quad_values, equations))

    # for scaling by curved geometric terms (not used by affine DGMultiMesh)
    flux_threaded = [[Trixi.allocate_nested_array(uEltype, nvars, (rd.Nq,), dg)
                      for _ in 1:NDIMS] for _ in 1:Threads.maxthreadid()]
    rotated_flux_threaded = [Trixi.allocate_nested_array(uEltype, nvars, (rd.Nq,), dg)
                             for _ in 1:Threads.maxthreadid()]

    # For backwards compatibility with older DGMulti code, the solution is included in the 
    # top-level cache. TODO: remove once DGMulti refactor is complete and stable.
    cache = (; md, weak_differentiation_matrices, lift_scalings, invJ, dxidxhatj,
             solution_container, u_values, u_face_values, flux_face_values,
             auxiliary_container, local_values_threaded, flux_threaded,
             rotated_flux_threaded)
    return cache
end

include("containers_manifold_covariant.jl")
include("dg_manifold_covariant.jl")
