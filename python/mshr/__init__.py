# -*- coding: utf-8 -*-
"""mshr package"""

# Copyright (C) 2017 Benjamin Kehlet

# This file is part of mshr.
#
# mshr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mshr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with mshr.  If not, see <http:#www.gnu.org/licenses/>.
#

import dolfin

from .cpp import Circle
from .cpp import Ellipse
from .cpp import Rectangle
from .cpp import Polygon
from .cpp import Sphere
from .cpp import Box
from .cpp import Cylinder
from .cpp import Cone
from .cpp import Tetrahedron
from .cpp import Surface3D
from .cpp import Ellipsoid
from .cpp import Rectangle
from .cpp import Extrude2D
from .cpp import CSGCGALDomain3D

from .cpp import CSGIntersection
from .cpp import CSGDifference
from .cpp import CSGUnion
from .cpp import CSGScaling
from .cpp import CSGTranslation
from .cpp import CSGRotation

from .cpp import CSGGeometries
from .cpp import UnitSphereMesh

from .cpp import CSGCGALMeshGenerator3D
from .cpp import TetgenMeshGenerator3D

from .cpp import CSGCGALDomain2D
from .cpp import CSGCGALMeshGenerator2D


from .cpp import generate_mesh
