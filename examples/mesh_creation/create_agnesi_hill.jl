using HOHQMesh

p = newProject("agnesi_hill", "src/meshes")
#all in km 
L = 409.6
H = 30.0
a = 16.0
h = 0.000016  

y_ll = (h * a^2) / (a^2 + (-L/2)^2) #B(X) from the paper needs to be shifted to the right by L/2
y_lr = (h * a^2) / (a^2 + (L/2)^2)

xEqn = "x(t) = t * $L"
yEqn = "y(t) = ($h * $a^2) / ($a^2 + (t*$L-($L/2))^2)"
zEqn = "z(t) = 0.0"

bottom = newParametricEquationCurve("bottom", xEqn, yEqn, zEqn)

xEqn = "x(t) = 0.0"
yEqn = "y(t) = $H + t * ($y_ll - $H)"
zEqn = "z(t) = 0.0"

left = newParametricEquationCurve("left", xEqn, yEqn, zEqn)


xEqn = "x(t) = $L"
yEqn = "y(t) = $y_lr + t * ($H - $y_lr)"
zEqn = "z(t) = 0.0"

right = newParametricEquationCurve("right", xEqn, yEqn, zEqn)

top = newEndPointsLineCurve("top", [L, H, 0.0], [0.0, H, 0.0])




addCurveToOuterBoundary!(p, bottom)
addCurveToOuterBoundary!(p, right)
addCurveToOuterBoundary!(p, top)
addCurveToOuterBoundary!(p, left)

setPolynomialOrder!(p, 2)
setPlotFileFormat!(p, "sem")

addBackgroundGrid!(p, [128.0, 300.0, 0.0]) #H1 in paper: dx=3.2 km, dz = 0.1 km
setMeshFileFormat!(p, "ABAQUS")
generate_mesh(p)