@muladd begin

#@inline function Trixi.get_node_coords


function init_auxiliary_node_variables!(auxiliary_variables, mesh::P4estMesh{2},
                                        equations::AbstractVariableCoefficientEquations{2,1}, dg,
                                        elements, velocity_field)
    (; node_coordinates) = elements
    (; aux_node_vars) = auxiliary_variables

    Trixi.@threaded for element in 1:Trixi.ncells(mesh)
        for j in eachnode(dg), i in eachnode(dg)
            # velocity field
            if !isnothing(velocity_field)
                x_node = Trixi.get_node_coords(node_coordinates, equations, dg, i, j,
                                                element)
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