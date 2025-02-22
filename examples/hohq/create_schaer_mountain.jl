using HOHQMesh

p = newProject("schaer_mountain_1000", "mesh")

y_top = 21000.0
x_ll = -25000.0
x_lr = 25000.0
width = x_lr - x_ll

h_c = 250.0
lambda_c = 4000.0
a_c = 5000.0

y_ll = h_c * exp(-((x_ll) / a_c)^2) * cos(pi * (x_ll) / lambda_c)^2 #height left side of the mountain
y_lr = h_c * exp(-((x_lr) / a_c)^2) * cos(pi * (x_lr) / lambda_c)^2 #height right side of the mountain



xEqn = "x(t) = $x_ll + t * $width"
yEqn = "y(t) = $h_c * exp(-(($x_ll + t * $width) / $a_c)^2) * cos(pi * ($x_ll + t * $width) / $lambda_c)^2"
zEqn = "z(t) = 0.0"

bottom = newParametricEquationCurve("bottom", xEqn, yEqn, zEqn)
#newEndPointsLineCurve("bottom", [-25000.0, 0.0,0.0], [25000.0,0.0 ,0.0])


#right   = newEndPointsLineCurve("right",[x_lr, y_lr, 0.0],  [x_lr, y_top, 0.0])
top     = newEndPointsLineCurve("top",  [x_lr, y_top, 0.0], [x_ll, y_top, 0.0])

xEqn = "x(t) = $x_lr"
yEqn = "y(t) = $y_lr + t * ($y_top - $y_lr)"
zEqn = "z(t) = 0.0"
right = newParametricEquationCurve("right", xEqn, yEqn, zEqn)

xEqn = "x(t) = $x_ll"
yEqn = "y(t) = $y_ll + (1 - t) * ($y_top - $y_ll)"
zEqn = "z(t) = 0.0"
left = newParametricEquationCurve("left", xEqn, yEqn, zEqn)

addCurveToOuterBoundary!(p, bottom)
addCurveToOuterBoundary!(p, right)
addCurveToOuterBoundary!(p, top)
addCurveToOuterBoundary!(p, left)

setPolynomialOrder!(p, 2)
setPlotFileFormat!(p, "sem")

addBackgroundGrid!(p, [1000.0, 840.0, 0.0])
setMeshFileFormat!(p, "ABAQUS")
generate_mesh(p)