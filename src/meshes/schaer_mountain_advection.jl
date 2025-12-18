using HOHQMesh

# ----------------------------------------------------------
# Project setup
# ----------------------------------------------------------
proj = newProject("schaer_mountain_advection","mesh")

setPolynomialOrder!(proj,1)
setMeshFileFormat!(proj,"ABAQUS")
setPlotFileFormat!(proj,"sem")

setName!(proj,"schaer_mountain_advection")

# Smoothing
setSmoothingStatus!(proj,"ON")
setSmoothingType!(proj,"LinearAndCrossbarSpring")

# Background grid
addBackgroundGrid!(proj,[50.0,25.0,0.0])

# ----------------------------------------------------------
# Refinement line (same as .control file)
# ----------------------------------------------------------
ref_line = newRefinementLine(
    "mountain_refinement",
    "smooth",
    [-18000.0,1500.0,0.0],
    [18000.0,1500.0,0.0],
    0.20,
    1500.0
)
addRefinementRegion!(proj, ref_line)

# ----------------------------------------------------------
# Outer boundary curves
# ----------------------------------------------------------

# ---- bottom_left: end-points line ----
bottom_left = newEndPointsLineCurve(
    "bottom_left",
    [-75000.0,0.0,0.0],
    [-25000.0,0.0,0.0]
)
addCurveToOuterBoundary!(proj, bottom_left)

# ---- bottom_left_connection: parametric line ----
bottom_left_con = newParametricEquationCurve(
    "bottom_left_connection",
    "x(t) = -25000.0",
    "y(t) = t * 9.600937856748462e-30",
    "z(t) = 0.0"
)
addCurveToOuterBoundary!(proj, bottom_left_con)

# ---- bottom: parametric mountain profile ----
bottom = newParametricEquationCurve(
    "bottom",
    "x(t) = -25000.0 + t*50000.0",
    """y(t) = cos(pi*(-25000.0 + t*50000.0)/8000.0)^2 *
              3000.0 *
              cos(pi*(-25000.0 + t*50000.0)/(2.0*25000.0))^2""",
    "z(t) = 0.0"
)
addCurveToOuterBoundary!(proj, bottom)

# ---- bottom_right_connection: parametric line ----
bottom_right_con = newParametricEquationCurve(
    "bottom_right_connection",
    "x(t) = 25000.0",
    "y(t) = 9.600937856748462e-30*(1 - t)",
    "z(t) = 0.0"
)
addCurveToOuterBoundary!(proj, bottom_right_con)

# ---- bottom_right: end-points line ----
bottom_right = newEndPointsLineCurve(
    "bottom_right",
    [25000.0,0.0,0.0],
    [75000.0,0.0,0.0]
)
addCurveToOuterBoundary!(proj, bottom_right)

# ---- right: end-points line ----
right = newEndPointsLineCurve(
    "right",
    [75000.0,0.0,0.0],
    [75000.0,25000.0,0.0]
)
addCurveToOuterBoundary!(proj, right)

# ---- top: end-points line ----
top = newEndPointsLineCurve(
    "top",
    [75000.0,25000.0,0.0],
    [-75000.0,25000.0,0.0]
)
addCurveToOuterBoundary!(proj, top)

# ---- left: end-points line ----
left = newEndPointsLineCurve(
    "left",
    [-75000.0,25000.0,0.0],
    [-75000.0,0.0,0.0]
)
addCurveToOuterBoundary!(proj, left)

# ----------------------------------------------------------
# Generate mesh
# ----------------------------------------------------------
generate_mesh(proj)
saveProject(proj)

println("Mesh generation complete!")
