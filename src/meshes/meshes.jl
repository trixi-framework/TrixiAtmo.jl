###############################################################################
# Mesh types for spherical geometry
###############################################################################

export face2global, face2cartesian, global2cartesian, global2face

"""
    Abstract type for a manifold (with boundary) embedded in Euclidean space
"""
abstract type AbstractManifold{dim, embedding_dim} end

export CubedSphereFace2D, cubed_sphere_connectivity
include("cubed_sphere_coupled_2d.jl")
