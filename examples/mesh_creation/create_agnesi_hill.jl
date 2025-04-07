using HOHQMesh

p = newProject("agnesi_hill", "src/meshes")

L =  409000.6
H = 30000.0
a = 16000.0
h = 0.016 

y_ll = (h * a^2) / (a^2 + (-(L/2))^2)
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

addBackgroundGrid!(p, [3200.0, 100.0, 0.0])
setMeshFileFormat!(p, "ABAQUS")
generate_mesh(p)