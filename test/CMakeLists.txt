# Copyright (C) 2014 Benjamin Kehlet
#
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
# along with mshr.  If not, see <http://www.gnu.org/licenses/>.

# This file defines regression tests. 

# Make sure python finds the mshr module.
# NOTE: Not sure this is the best way. Is saves the environment at "cmake
# time" and applies it when running the test. Changing PYTHONPATH
# later will not have effect. Not sure how to add tp PYTHONPATH at
# runtime.
set(PYTHON_ENVIR "$ENV{PYTHONPATH};DOLFIN_NOPLOT=True")


############# A python program that just imports dolfin and mshr ##############################
add_test("Python-DummyImport" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/dummy.py")
set_property(TEST "Python-DummyImport" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
)

############# A union which is particularly prone to roundoff errors ##########################
add_test("Python-DegenerateRemoval" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/degenerate_removal.py")
set_property(TEST "Python-DegenerateRemoval" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
)

############# Test the csg operators ##########################################################
add_test("Python-CSGOperators" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/test-csg.py")
set_property(TEST "Python-CSGOperators" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
)

############# Test the fuzzy point map ########################################################
# This test is disabled since the Surface3D class is about to be removed.
#add_test("Python-FuzzyPointMap" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/test-fuzzypointmap.py")
#set_property(TEST "Python-FuzzyPointMap" PROPERTY
#    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
#)

############# Test the mesh generation ########################################################
add_test("Python-MeshGeneration" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/test-mesh-generation.py")
set_property(TEST "Python-MeshGeneration" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
)

############# Test the predefined meshes ##########################################################
add_test("Python-Meshes" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/test-meshes.py")
set_property(TEST "Python-Meshes" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
)
############# Test the ASCFileReader ##########################################################
add_test("Python-ASCFileReader" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/test-ASCFileReader.py")
set_property(TEST "Python-ASCFileReader" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
)

############# Test the ASCFileReader ##########################################################
add_test("Python-NumSegments2D" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/test-num-segments-2d.py")
set_property(TEST "Python-NumSegments2D" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
)

############# Test the CSG predicates (2D) ##########################################################
add_test("Python-CSGPredicates2D" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/test-csg-predicates.py")
set_property(TEST "Python-CSGPredicates2D" PROPERTY
  ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
  )

############# Test CSG primitves (2D) ##########################################################
add_test("Python-CSGPrimitives2D" "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/test/test-csg-primitives-2d.py")
set_property(TEST "Python-CSGPredicates2D" PROPERTY
  ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
)

############# Try meshing some surface files found in the cgal source ########################
# TODO: Print mesh quality measure from program and check the value
set(CGAL_DATA_DIR "${CMAKE_SOURCE_DIR}/3rdparty/CGAL/demo/Polyhedron/data")

# Files that should pass
# Note elephant.off is valid a should passes but takes a very long (so test slave may time out)
set(TESTFILES_VALID "anchor.off;couplingdown.off;cross.off;cube.off;dragknob.off;ellipsoid.off;handle.off;icosahedron.off;joint.off;knot1.off;knot2.off;pinion.off;pipe.off;rotor.off;sphere.off;spool.off;star.off;translated-cube.off;tripod.off")
set(TESTFILES_SELFINTERSECTING "bones.off;cow.off;man.off;oblong.ogg")
set(TESTFILES_NOTCLOSED "cube-ouvert.off;mushroom.off")
set(TESTFILES_NOTTRIANGULAR "pyramid.off")

foreach(CURRENT_DATA_FILE IN LISTS TESTFILES_VALID)
  message(STATUS "Adding test: ${CURRENT_DATA_FILE}")
  add_test(FilePass-cgal-${CURRENT_DATA_FILE} ../mshrable -s -b cgal --check-mesh ${CGAL_DATA_DIR}/${CURRENT_DATA_FILE})
  add_test(FilePass-tetgen-${CURRENT_DATA_FILE} ../mshrable -s -b tetgen --check-mesh ${CGAL_DATA_DIR}/${CURRENT_DATA_FILE})
endforeach(CURRENT_DATA_FILE)


##############  CPP tests
add_executable(cppTestCSGGeometries test-csggeometries.cpp)
target_link_libraries(cppTestCSGGeometries mshr)
add_test(CPP-CSGGeometries cppTestCSGGeometries)


############ Run the demos as regression test ################################################
file(GLOB PYTHONDEMOS "${CMAKE_SOURCE_DIR}/demo/python/*.py")
foreach(CURRENT_DEMO IN LISTS PYTHONDEMOS)
  add_test("PythonDemo-${CURRENT_DEMO}" "${PYTHON_EXECUTABLE}" "${CURRENT_DEMO}")
  set_property(TEST "PythonDemo-${CURRENT_DEMO}" PROPERTY
    ENVIRONMENT "PYTHONPATH=${PYTHON_ENVIR}"
  )

endforeach()
