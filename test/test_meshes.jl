@testitem "Meshes: P4estMesh connectivities are complete" setup=[Setup] tags=[:meshes] begin
    # Without all corners (and edges in 3D) connecting trees, a single call to
    # `p4est_balance` does not always produce a balanced forest, see
    # https://github.com/cburstedde/p4est/issues/112.
    # Refining the quadrant at the origin of each tree recursively to level 5 is enough to
    # trigger this for incomplete connectivities.
    function refine_corner_2d(p4est, which_tree, quadrant)
        q = unsafe_load(quadrant)
        return Cint(q.x == 0 && q.y == 0 && q.level < 5)
    end
    refine_corner_2d_c = @cfunction(refine_corner_2d, Cint,
                                    (Ptr{Trixi.p4est_t}, Trixi.p4est_topidx_t,
                                     Ptr{Trixi.p4est_quadrant_t}))

    function refine_corner_3d(p8est, which_tree, quadrant)
        q = unsafe_load(quadrant)
        return Cint(q.x == 0 && q.y == 0 && q.z == 0 && q.level < 5)
    end
    refine_corner_3d_c = @cfunction(refine_corner_3d, Cint,
                                    (Ptr{Trixi.p8est_t}, Trixi.p4est_topidx_t,
                                     Ptr{Trixi.p8est_quadrant_t}))

    # The numbers of vertices, edges, and corners follow from the topology of the
    # meshes, e.g., the Euler characteristic of the sphere for the 2D meshes. Corners
    # (edges) are only stored if they connect trees that are not already connected via
    # faces (and edges), e.g., not at the vertices of the cube in the 2D cubed sphere.
    meshes_2d = [
        (P4estMeshCubedSphere2D(1, 1.0; polydeg = 1), 8, 0),
        (P4estMeshCubedSphere2D(2, 1.0; polydeg = 1), 26, 18),
        (P4estMeshQuadIcosahedron2D(2, 1.0; polydeg = 1), 242, 222)
    ]
    for (mesh, n_vertices, n_corners) in meshes_2d
        connectivity = mesh.p4est.connectivity
        @test connectivity.num_vertices[] == n_vertices
        @test connectivity.num_corners[] == n_corners

        # Call `p4est_balance` directly since `Trixi.balance!` may call it twice
        Trixi.refine_p4est!(mesh.p4est, true, refine_corner_2d_c, C_NULL)
        Trixi.p4est_balance(mesh.p4est, Trixi.P4EST_CONNECT_FACE, C_NULL)
        @test Trixi.p4est_is_balanced(mesh.p4est, Trixi.P4EST_CONNECT_FACE) == 1
    end

    for keep_columns_together in (false, true)
        mesh = P4estMeshCubedSphereTopography(2, 2, 1.0, 0.5; polydeg = 1,
                                              keep_columns_together)
        connectivity = mesh.p4est.connectivity
        @test connectivity.num_vertices[] == 78
        @test connectivity.num_edges[] == 84
        @test connectivity.num_corners[] == 18

        Trixi.refine_p4est!(mesh.p4est, true, refine_corner_3d_c, C_NULL)
        Trixi.p8est_balance(mesh.p4est, Trixi.P8EST_CONNECT_FACE, C_NULL)
        @test Trixi.p8est_is_balanced(mesh.p4est, Trixi.P8EST_CONNECT_FACE) == 1
    end
end

@testitem "Meshes: complete_connectivity keeps face connections" setup=[Setup] tags=[:meshes] begin
    # `complete_connectivity` throws an error if the face connections change, so it
    # suffices to complete the face-only connectivities. Here, we additionally check that
    # the completed connectivity contains all neighbors of the 2x2 corner of a brick.
    connectivity_brick = Trixi.p4est_connectivity_new_brick(2, 2, 0, 0)
    connectivity = TrixiAtmo.complete_connectivity(connectivity_brick)
    connectivity_pw = Trixi.PointerWrapper(connectivity)
    @test connectivity_pw.num_vertices[] == 9
    @test connectivity_pw.num_corners[] == 1
    @test unsafe_wrap(Array, connectivity_pw.corner_to_tree, 4) == [0, 1, 2, 3]
    @test unsafe_wrap(Array, connectivity_pw.corner_to_corner, 4) == [3, 2, 1, 0]
    Trixi.p4est_connectivity_destroy(connectivity)

    # A single tree that is periodic in one direction touches itself
    connectivity = Trixi.p4est_connectivity_new_brick(1, 1, 1, 0)
    @test_throws ErrorException TrixiAtmo.complete_connectivity(connectivity)
end
