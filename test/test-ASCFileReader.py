from dolfin import *
from mshr import *

import os
import tempfile

def save_surface(filename) :
        geometry = Sphere(Point(0.0,0.0,0.0),10.)
        domain = CSGCGALDomain3D(geometry)
        domain.save(filename)

def load_surface(filename) :
        surf  = Surface3D(filename)
        domain = CSGCGALDomain3D(surf)
        mesh = generate_mesh(domain, 10)

fd, temp_path = tempfile.mkstemp(suffix='.asc')

save_surface(temp_path)
load_surface(temp_path)

os.close(fd)
os.remove(temp_path)
