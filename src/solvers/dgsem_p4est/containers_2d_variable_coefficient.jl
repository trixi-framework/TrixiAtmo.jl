@muladd begin

    #@inline function Trixi.get_node_coords


    function init_auxiliary_node_variables!(
        auxiliary_variables,
        mesh::P4estMesh{2},
        equations::AbstractVariableCoefficientEquations{2,1},
        dg,
        elements,
        velocity_field,
    )
        (; node_coordinates) = elements
        (; aux_node_vars) = auxiliary_variables

        Trixi.@threaded for element = 1:Trixi.ncells(mesh)
            for j in eachnode(dg), i in eachnode(dg)
                # velocity field
                if !isnothing(velocity_field)
                    x_node = Trixi.get_node_coords(
                        node_coordinates,
                        equations,
                        dg,
                        i,
                        j,
                        element,
                    )
                    aux_node_vars[1, i, j, element] = velocity_field(x_node)[1]
                    aux_node_vars[2, i, j, element] = velocity_field(x_node)[2]
                else
                    aux_node_vars[1, i, j, element] = zero(eltype(aux_node_vars))
                    aux_node_vars[2, i, j, element] = zero(eltype(aux_node_vars))
                end
            end
        end
        return nothing
    end

    function Trixi.init_elements(
        mesh::P4estMesh{2,2,RealT},
        equations::VariableCoefficientAdvectionEquation2D,
        basis,
        metric_terms,
        ::Type{uEltype},
    ) where {RealT<:Real,uEltype<:Real}
        nelements = Trixi.ncells(mesh)

        _node_coordinates = Vector{RealT}(undef, 2 * nnodes(basis)^2 * nelements)
        node_coordinates = unsafe_wrap(
            Array,
            pointer(_node_coordinates),
            (2, ntuple(_ -> nnodes(basis), 2)..., nelements),
        )

        _jacobian_matrix = Vector{RealT}(undef, 2^2 * nnodes(basis)^2 * nelements)
        jacobian_matrix = unsafe_wrap(
            Array,
            pointer(_jacobian_matrix),
            (2, 2, ntuple(_ -> nnodes(basis), 2)..., nelements),
        )

        _contravariant_vectors = similar(_jacobian_matrix)
        contravariant_vectors =
            unsafe_wrap(Array, pointer(_contravariant_vectors), size(jacobian_matrix))

        _inverse_jacobian = Vector{RealT}(undef, nnodes(basis)^2 * nelements)
        inverse_jacobian = unsafe_wrap(
            Array,
            pointer(_inverse_jacobian),
            (ntuple(_ -> nnodes(basis), 2)..., nelements),
        )

        _surface_flux_values = Vector{uEltype}(
            undef,
            nvariables(equations) * nnodes(basis)^(2 - 1) * (2 * 2) * nelements,
        )
        surface_flux_values = unsafe_wrap(
            Array,
            pointer(_surface_flux_values),
            (
                nvariables(equations),
                ntuple(_ -> nnodes(basis), 2 - 1)...,
                2 * 2,
                nelements,
            ),
        )

        elements = Trixi.P4estElementContainer{2,RealT,uEltype,2 + 1,2 + 2,2 + 3}(
            node_coordinates,
            jacobian_matrix,
            contravariant_vectors,
            inverse_jacobian,
            surface_flux_values,
            _node_coordinates,
            _jacobian_matrix,
            _contravariant_vectors,
            _inverse_jacobian,
            _surface_flux_values,
        )

        Trixi.init_elements!(elements, mesh, basis)
        return elements
    end

end #muladd
