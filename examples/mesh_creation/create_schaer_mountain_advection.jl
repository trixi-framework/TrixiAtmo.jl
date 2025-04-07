using HOHQMesh

p = newProject("schaer_mountain_advection", "src/meshes")

y_top = 25000.0
x_ll = -75000.0
x_lr = 75000.0

h_0 = 3000.0
lambda = 8000.0
a = 25000.0
a_neg = -25000.0
width = a - a_neg

y_neg_a = cos(pi * a_neg / lambda)^2.0 * h_0 * cos(pi * a_neg / (2.0 * a))^2.0
y_a = cos(pi * a / lambda)^2.0 * h_0 * cos(pi * a / (2.0 * a))^2.0



xEqn = "x(t) = $a_neg + t * $width"
yEqn = "y(t) = cos(pi*($a_neg+t*$width)/$lambda)^2.0*$h_0*cos(pi*($a_neg+t*$width)/(2.0*$a))^2.0"
zEqn = "z(t) = 0.0"

bottom = newParametricEquationCurve("bottom", xEqn, yEqn, zEqn)


xEqn = "x(t) = $a_neg"
yEqn = "y(t) = t * $y_neg_a"
zEqn = "z(t) = 0.0"

bottom_left_connection = newParametricEquationCurve("bottom_left_connection", xEqn, yEqn, zEqn)


xEqn = "x(t) = $a"
yEqn = "y(t) = $y_a - t * $y_a"
zEqn = "z(t) = 0.0"

bottom_right_connection = newParametricEquationCurve("bottom_right_connection", xEqn, yEqn, zEqn)

bottom_left = newEndPointsLineCurve("bottom_left", [x_ll, 0.0, 0.0], [a_neg, 0.0, 0.0])
bottom_right = newEndPointsLineCurve("bottom_right", [a, 0.0, 0.0], [x_lr, 0.0, 0.0])


#right   = newEndPointsLineCurve("right",[x_lr, y_lr, 0.0],  [x_lr, y_top, 0.0])
top = newEndPointsLineCurve("top", [x_lr, y_top, 0.0], [x_ll, y_top, 0.0])


left = newEndPointsLineCurve("left", [x_ll, y_top, 0.0], [x_ll, 0.0, 0.0])
right = newEndPointsLineCurve("right", [x_lr, 0.0, 0.0], [x_lr, y_top, 0.0])


addCurveToOuterBoundary!(p, bottom_left)
addCurveToOuterBoundary!(p, bottom_left_connection)
addCurveToOuterBoundary!(p, bottom)
addCurveToOuterBoundary!(p, bottom_right_connection)
addCurveToOuterBoundary!(p, bottom_right)
addCurveToOuterBoundary!(p, right)
addCurveToOuterBoundary!(p, top)
addCurveToOuterBoundary!(p, left)

setPolynomialOrder!(p, 2)
setPlotFileFormat!(p, "sem")

addBackgroundGrid!(p, [1000.0, 840.0, 0.0])
setMeshFileFormat!(p, "ABAQUS")
generate_mesh(p)