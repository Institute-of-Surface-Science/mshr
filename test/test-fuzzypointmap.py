import mshr
import tempfile, os, sys

# A simple unit cube but with some vertices perturbed slightly to test the FuzzyPointMap
cube = """
solid ascii
facet normal 0 0 0
  outer loop
    vertex 1e-12 0 1
    vertex 0 1 0
    vertex 0 -1e-11 0
  endloop
endfacet

facet normal 0 0 0
outer loop
  vertex 0 0 1
  vertex 0 1 1
  vertex 0 1 1e-13
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 0 0 1
  vertex 1 0 1
  vertex 0 1 1
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 1 0 1
  vertex 1 1 1
  vertex 0 1 1
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 1 0 1
  vertex 1 0 0
  vertex 1 1 1
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 1 0 0
  vertex 1 1 0
  vertex 1 1 1
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 1 0 0
  vertex 0 0 0
  vertex 1 1 0
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 0 0 0
  vertex 0 1 0
  vertex 1 1 0
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 1 1 1
  vertex 1 1 0
  vertex 0 1 1
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 1 1 0
  vertex 0 1 0
  vertex 0 1 1
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 0 0 1
  vertex 0 0 0
  vertex 1 0 1
endloop
endfacet
facet normal 0 0 0
outer loop
  vertex 0 0 0
  vertex 1 0 0
  vertex 1 0 1
endloop
endfacet
endsolid
"""

# In python 3 convert the text from unicode to a "classic" str
if sys.version_info[0] >= 3 :
    cube = cube.encode("ascii")

fd, filename = tempfile.mkstemp(suffix=".stl")
os.write(fd, cube)
os.close(fd)

s = mshr.Surface3D(filename)
d = mshr.CSGCGALDomain3D(s)

os.remove(filename)

assert d.num_holes() == 0
