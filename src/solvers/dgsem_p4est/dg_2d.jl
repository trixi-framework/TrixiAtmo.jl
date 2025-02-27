@muladd begin
    #! format: noindent

@inline function Trixi.weak_form_kernel!



end

@inline function Trixi.calc_interface_flux!(surface_flux_values, mesh::P4estMesh{2},
    nonconservative_terms::False,
    equations::AbstractCovariantEquations{2},
    surface_integral, dg::DG, cache,
    interface_index, normal_direction,
    primary_node_index,
    primary_direction_index,
    primary_element_index,
    secondary_node_index,
    secondary_direction_index,
    secondary_element_index)



end

@inline function Trixi.calc_interface_flux!(surface_flux_values, mesh::P4estMesh{2},
    nonconservative_terms::True,
    equations::AbstractCovariantEquations{2},
    surface_integral, dg::DG, cache,
    interface_index, normal_direction,
    primary_node_index,
    primary_direction_index,
    primary_element_index,
    secondary_node_index,
    secondary_direction_index,
    secondary_element_index)


end