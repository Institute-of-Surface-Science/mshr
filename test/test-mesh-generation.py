import mshr
from dolfin import *

# issue 37
g = mshr.Rectangle(Point(0.0, 0.0), Point(2.2, .41)) - mshr.Circle(Point(.2, .2), .05, 40)
m = mshr.generate_mesh(g, 50)

# issue 41 (failed only in parallel)
c = mshr.Extrude2D(mshr.Circle(Point(0, 0, 0), 1), 1)
m = mshr.generate_mesh(c, 10)
