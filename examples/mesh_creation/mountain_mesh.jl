using HOHQMesh
#using GLMakie

project = newProject("Mountain", "src/meshes")  # project name, output folder

mountainHeight = 400      # height of montain 400 m
mountainHalfWidth = 1000  # half width of mountain 1 km
x_ll = -20000.0           # domain width 40 km
x_lr = 20000.0
width = x_lr - x_ll
y_top = 16000.0           # domain height 16 km

# lower left and right values according to parametrization
y_ll = mountainHeight / ( 1 + (x_ll / mountainHalfWidth)^2 )
y_lr = mountainHeight / ( 1 + (x_lr / mountainHalfWidth)^2 )
#y_ll = 0.0
#y_lr = 0.0
# println( "y_bottom = ", y_ll, ", ", y_lr)

# straight line for bottom
# bottom  = newEndPointsLineCurve("bottom",  [x_ll, y_ll,  0.0], [x_lr, y_lr,  0.0])

# parametrized bottom curve, t ∈ [0,1]
xEqn = "x(t) = $x_ll + t * $width"
yEqn = "y(t) = $mountainHeight / ( 1 + (($x_ll + t * $width) / $mountainHalfWidth)^2 )"
zEqn = "z(t) = 0.0"
bottom = newParametricEquationCurve("bottom", xEqn, yEqn, zEqn)

# straight lines for two boundaries
right   = newEndPointsLineCurve("right",   [x_lr, y_lr,  0.0], [x_lr, y_top, 0.0])
top     = newEndPointsLineCurve("top",     [x_lr, y_top, 0.0], [x_ll, y_top, 0.0])

# workaround to use a parametric equation for a straight line. Has slightly different round-off
# errors that are not as sensitive
xEqn = "x(t) = $x_ll"
yEqn = "y(t) = $y_ll + (1 - t) * ($y_top - $y_ll)"
zEqn = "z(t) = 0.0"
left = newParametricEquationCurve("left", xEqn, yEqn, zEqn)
# left = newEndPointsLineCurve("left",    [x_ll, y_top, 0.0], [x_ll, y_ll,  0.0])
x
# use the four curves as outer boundary
for curve in [right, top, left, bottom]
    addCurveToOuterBoundary!(project, curve)
end

# set the desired resolution in x and y direction
addBackgroundGrid!(project, [2000.0, 2000.0, 0.0])

# have a look
#plotProject!(project, MODEL+GRID)

# generate
#@info "Press enter to generate the mesh and update the plot."
#readline()
generate_mesh(project)

@info "Press enter to quit."
readline()
