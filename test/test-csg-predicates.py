from dolfin import *
import mshr
import math


r = mshr.Rectangle(Point(0,-.5), Point(4,.5))

# Trivial inside outside
assert r.inside(Point( .5,  .0))
assert not r.inside(Point(1.5, 1.5))

rotated = mshr.CSGRotation(r, math.pi/2)
assert not rotated.inside(Point(1.5, 1.5))
assert rotated.inside(Point(0., 1.5))
assert not rotated.inside(Point(0., 4.5))
assert not rotated.inside(Point(0., -.5))

# Translation
translated = mshr.CSGTranslation(rotated, Point(2, 1))
assert translated.inside(Point(2, 1.5))
assert not translated.inside(Point(0, 1.5))

# Circle
c = mshr.Circle(Point(2, 1), .5)
assert c.inside(Point(2, 1))
assert not c.inside(Point(3, 1))

c_rotated = mshr.CSGRotation(c, Point(2,1), math.pi)
assert c_rotated.inside(Point(2,1))
assert not c_rotated.inside(Point(3, 1))

c_rotated_2 = mshr.CSGRotation(c_rotated, pi/2)
assert not c_rotated_2.inside(Point(2,1))
assert c_rotated_2.inside(Point(-1,2))

# Rectangle and box with the vertices not ordered
assert mshr.Rectangle(Point(0,0),Point(-1,-1)).inside(Point(-0.5,-0.5))
#assert mshr.Box(Point(0,0,0),Point(-1,-1,-1)).inside(Point(-0.5,-0.5,-0.5))
