from dolfin import *
import mshr

# The ellipse is completely contained in the circle and will not be
# visible in the resulting mesh, but it is so small that with
# mesh_resolution = 10 that when converting it to a polygon the number
# of segments in the polygon will be 0 (ie. the polygon is degenerate)
# if not handled correctly.
c = mshr.Circle(Point(0, 0), 1.5)
e = mshr.Ellipse(Point(.05, 0), .01, .05)

mesh = mshr.generate_mesh(c+e, 10)

