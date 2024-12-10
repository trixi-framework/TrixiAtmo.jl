using HOHQMesh

p = newProject("schaer_mountain_boundary", "out")

bottom = newParametricEquationCurve("bottom", "x(t)=t*50000-25000", "y(t)=250.0 * exp(-((t*50000-25000) /  5000.0)^2) * cos(pi * (t*50000-25000) / 4000.0)^2", "z(t)=0.0") #newEndPointsLineCurve("bottom", [-25000.0, 0.0,0.0], [25000.0,0.0 ,0.0])

h=250.0 * exp(-((1*50000-25000) /  5000.0)^2) * cos(pi * (1*50000-25000) / 4000.0)^2

right = newEndPointsLineCurve("right", [25000.0,h,0.0], [25000.0,21000.0 ,0.0])
top = newEndPointsLineCurve("top", [25000.0,21000.0,0.0], [-25000.0,21000.0,0.0])
left = newEndPointsLineCurve("left", [-25000.0,21000.0,0.0], [-25000,h,0.0])

addCurveToOuterBoundary!(p, bottom)
addCurveToOuterBoundary!(p, right)
addCurveToOuterBoundary!(p, top)
addCurveToOuterBoundary!(p, left) 

setPolynomialOrder!(p, 2)
setPlotFileFormat!(p, "sem")

addBackgroundGrid!(p, [1000.0, 840.0, 0.0])

generate_mesh(p)