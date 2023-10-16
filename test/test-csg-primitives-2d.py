from dolfin import *
from mshr import *

# Polygon with list of points as argument
p = Polygon([Point(0, 0), Point(1, 1), Point(0, 1)])

# Polygon with tuple of points as argument
p = Polygon( (Point(0, 0), Point(1, 1), Point(0, 1)) )

