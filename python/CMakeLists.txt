cmake_minimum_required(VERSION 3.5.0)


# This script is used solely for getting the paths to the dependencies
# of the bindings (ie. pybind11, dolfin and mshr).
# This is written to a json file which is read from setup.py.

PROJECT(mshr_pybind11_config)


#   pybind11_FOUND - true if pybind11 and all required components found on the system
#   pybind11_VERSION - pybind11 version in format Major.Minor.Release
#   pybind11_INCLUDE_DIRS - Directories where pybind11 and python headers are located.
#   pybind11_INCLUDE_DIR - Directory where pybind11 headers are located.
#   pybind11_DEFINITIONS - Definitions necessary to use pybind11, namely USING_pybind11.
#   pybind11_LIBRARIES - compile flags and python libraries (as needed) to link against.
#   pybind11_LIBRARY - empty.
#   CMAKE_MODULE_PATH - appends location of accompanying FindPythonLibsNew.cmake and
#                       pybind11Tools.cmake modules.
find_package(pybind11 CONFIG HINTS ${PYBIND11_DIR} ${PYBIND11_ROOT}
  $ENV{PYBIND11_DIR} $ENV{PYBIND11_ROOT})


find_package(DOLFIN)

find_package(mshr)

configure_file("config.json.in" "config.json")


