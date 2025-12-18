using HOHQMesh

proj = newProject("schaer_mountain_advection3","mesh")

setPolynomialOrder!(proj,1)
setMeshFileFormat!(proj,"ABAQUS")
setPlotFileFormat!(proj,"sem")
setName!(proj,"schaer_mountain_advection3")

setSmoothingStatus!(proj,"ON")
setSmoothingType!(proj,"LinearAndCrossbarSpring")
setSmoothingIterations!(proj,2)


addBackgroundGrid!(proj,[500.0,250.0,0.0])

ref_line = newRefinementLine("mountain_refinement","smooth",[-18000.0,1500.0,0.0],[18000.0,1500.0,0.0],0.20,1500.0)
addRefinementRegion!(proj,ref_line)
x_eqn_bottom = "x(t) = -75000.0 + t * 150000.0"
y_eqn_bottom = "y(t) = 3000.0 * exp(-((-75000.0 + t * 150000.0) / 9511.0)^2) * cos(pi * (-75000.0 + t * 150000.0) / 8000.0)^2"
z_eqn_bottom = "z(t) = 0.0"


bottom = newParametricEquationCurve("bottom",x_eqn_bottom,y_eqn_bottom,z_eqn_bottom)
addCurveToOuterBoundary!(proj,bottom)

right = newParametricEquationCurve("right","x(t)=75000.0","y(t)=4.336711471172962322e-25 + t*(15000.0 - 4.336711471172962322e-25)","z(t)=0.0")
addCurveToOuterBoundary!(proj,right)

top = newEndPointsLineCurve("top",[75000.0,15000.0,0.0],[-75000.0,15000.0,0.0])
addCurveToOuterBoundary!(proj,top)

left = newParametricEquationCurve("left","x(t)=-75000.0","y(t)=4.336711471172962322e-25 + (1 - t)*(15000.0 - 4.336711471172962322e-25)","z(t)=0.0")

addCurveToOuterBoundary!(proj,left)
generate_mesh(proj)
saveProject(proj)
println("Mesh generation complete!")
