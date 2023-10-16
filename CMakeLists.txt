# Require CMake 3.5
cmake_minimum_required(VERSION 3.5)

project( MSHR )
set(MSHR_VERSION_RELEASE 1)
set(MSHR_VERSION_MAJOR "2019")
set(MSHR_VERSION_MINOR "1")
set(MSHR_VERSION_MICRO "0")
set(MSHR_VERSION "${MSHR_VERSION_MAJOR}.${MSHR_VERSION_MINOR}.${MSHR_VERSION_MICRO}")
if (NOT MSHR_VERSION_RELEASE)
  set(MSHR_VERSION "${MSHR_VERSION}.dev0")
endif()

# CGAL setup
option(USE_SYSTEM_CGAL "Do not build CGAL, but use an existing build instead." OFF)
if (USE_SYSTEM_CGAL)
  find_package(CGAL 4.12 CONFIG REQUIRED)
endif()

# Borrow some cmake modules from cgal
if (USE_SYSTEM_CGAL)
  set(CMAKE_MODULE_PATH "${CGAL_MODULES_DIR}")
else()
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/3rdparty/CGAL/cmake/modules")
endif()

# Add cmake directory to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(USE_SYSTEM_TETGEN "Do not build tetgen, but use an existing build instead." OFF)
if (USE_SYSTEM_TETGEN)
  find_package(TetGen)
endif()

# Helper macro for testing if particular value is contained in list
# Taken from http://www.cmake.org/Wiki/CMakeMacroListOperations#LIST_CONTAINS
MACRO(LIST_CONTAINS var value)
  SET(${var})
  FOREACH (value2 ${ARGN})
    IF (${value} STREQUAL ${value2})
      SET(${var} TRUE)
    ENDIF (${value} STREQUAL ${value2})
  ENDFOREACH (value2)
ENDMACRO(LIST_CONTAINS)

# Use C++11
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall")

# Boost
# This is workaround to avoid that find_package(Boost)
# picks up th ewrong boost when a hint is given
set(BOOST_ROOT $ENV{BOOST_DIR} $ENV{BOOST_HOME})
if (BOOST_ROOT)
  set(Boost_NO_SYSTEM_PATHS on)
endif()

# Prevent FindBoost.cmake from looking for system Boost{foo}.cmake files
set(Boost_NO_BOOST_CMAKE true)
find_package( Boost REQUIRED system filesystem program_options )
# Save this because it will be overwritten by find_package(DOLFIN)
set(Boost_MSHR_LIBRARIES ${Boost_LIBRARIES})
include_directories("${Boost_INCLUDE_DIRS}")

#  GMP_INCLUDE_DIR       - the GMP include directory
#  GMP_LIBRARIES_DIR     - directory where the GMP libraries are located
#  GMP_LIBRARIES         - Link these to use GMP
find_package(GMP REQUIRED)
include_directories("${GMP_INCLUDE_DIR}")

# Try to find the MPFR libraries
# MPFR_FOUND - system has MPFR lib
# MPFR_INCLUDE_DIR - the MPFR include directory
# MPFR_LIBRARIES_DIR - Directory where the MPFR libraries are located
# MPFR_LIBRARIES - the MPFR libraries
# MPFR_IN_CGAL_AUXILIARY - TRUE if the MPFR found is the one distributed with CGAL in the auxiliary folder
find_package(MPFR REQUIRED)
include_directories("${MPFR_INCLUDE_DIR}")

find_package(Eigen3 REQUIRED)
include_directories("${EIGEN3_INCLUDE_DIR}")

# use CGAL in header-only mode
option (CGAL_HEADER_ONLY "Use CGAL in header-only mode" ON)
if (CGAL_HEADER_ONLY)
  add_definitions( -DCGAL_HEADER_ONLY )
endif()

add_subdirectory(3rdparty)
#include_directories(BEFORE ${CGAL_INCLUDE_DIR})
include_directories(BEFORE ${EXTERNAL_INCLUDE_DIRS})
add_definitions("${EXTERNAL_DEFINITIONS}")

find_package(DOLFIN REQUIRED)
include_directories(${DOLFIN_PYTHON_INCLUDE_DIRS})
include(${DOLFIN_USE_FILE})

# Set installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(SHARE_DIR "share/mshr" CACHE PATH "Shared data installation directory.")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()


# include for local header directory
include_directories( BEFORE include )

file( GLOB_RECURSE SOURCES src/*.cpp )

add_library(mshr SHARED ${SOURCES})
if (NOT USE_SYSTEM_CGAL)
  add_dependencies( mshr CGAL )
endif()

# Append the library version information to the library target properties
option(MSHR_WITH_LIBRARY_VERSION "Build with library version information." ON)
if (MSHR_WITH_LIBRARY_VERSION)
  string(REPLACE "+" "" MSHR_LIBRARY_VERSION ${MSHR_VERSION})
  # This setting of SOVERSION assumes that any API change
  # will increment either the minor or major version number.
  set(MSHR_LIBRARY_PROPERTIES ${MSHR_LIBRARY_PROPERTIES}
    VERSION ${MSHR_LIBRARY_VERSION}
    SOVERSION ${MSHR_VERSION_MAJOR}.${MSHR_VERSION_MINOR}
  )
endif()
set_target_properties(mshr PROPERTIES ${MSHR_LIBRARY_PROPERTIES})

# Link the executable to (static) CGAL libraries
target_link_libraries( mshr ${Boost_MSHR_LIBRARIES}
                            ${GMP_LIBRARIES}
                            ${MPFR_LIBRARIES}
			    tet
                            ${DOLFIN_LIBRARIES}
                     )

# Fix CGAL's random generator (usefull to reproduce bugs)
option(INIT_CGAL_RANDOM "Use fixed seed for CGAL's pseudo random generator" OFF)
if (INIT_CGAL_RANDOM)
  add_definitions( -DINIT_RANDOM_GENERATOR=0 )
endif()

#
option(ENABLE_EXPERIMENTAL "Enable experimental code" OFF)
if(ENABLE_EXPERIMENTAL)
  message(STATUS "Enabling experimental code")
  add_definitions( -DMSHR_ENABLE_EXPERIMENTAL )
endif()


export(PACKAGE mshr)

target_link_libraries( mshr
  ${DOLFIN_LIBRARIES}
  ${Boost_MSHR_LIBRARIES}
  )

# install library
install(TARGETS mshr
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  )

option(ENABLE_MSHRABLE "Enable building the command line client, mshrable" OFF)
if (ENABLE_MSHRABLE)
  add_executable( mshrable app/mshrable.cpp )

  target_link_libraries( mshrable
    mshr
    ${DOLFIN_LIBRARIES}
    ${Boost_MSHR_LIBRARIES}
    )

  # install app
  install(TARGETS mshrable
    RUNTIME DESTINATION ${INSTALL_BIN_DIR}
    LIBRARY DESTINATION ${INSTALL_LIB_DIR}
    ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
    )
else()
    message(STATUS "Disabling building command line client mshrable")
endif()

# install header files
install(DIRECTORY include/ DESTINATION ${INSTALL_INCLUDE_DIR})

file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")

### add cmake configuration file ###

# common for both config files
set(CONF_EXTERNAL_INCLUDE_DIRS "${DOLFIN_PYTHON_INCLUDE_DIRS};${DOLFIN_INCLUDE_DIRS};${DOLFIN_3RD_PARTY_INCLUDE_DIRS}")
set(CONF_EXTERNAL_LIBRARIES "${DOLFIN_PYTHON_LIBRARIES};${DOLFIN_LIBRARIES};${DOLFIN_3RD_PARTY_LIBRARIES}")
set(CONF_CXX_DEFINITIONS "${DOLFIN_CXX_DEFINITIONS}")
string(REPLACE "\"" "\\\"" CONF_CXX_DEFINITIONS "${CONF_CXX_DEFINITIONS}")
set(CONF_CXX_FLAGS "${DOLFIN_CXX_FLAGS}")

# for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")
set(CONF_LIBRARIES_DIRS "${PROJECT_BINARY_DIR}")

configure_file(mshrConfig.cmake.in "${PROJECT_BINARY_DIR}/mshrConfig.cmake" @ONLY)
configure_file(mshr-config.cmake.in "${PROJECT_BINARY_DIR}/mshr-config.cmake" @ONLY)

# ... and for the install tree
set(CONF_INCLUDE_DIRS "${INSTALL_INCLUDE_DIR}")
set(CONF_LIBRARIES_DIRS "${INSTALL_LIB_DIR}")
configure_file(mshrConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mshrConfig.cmake" @ONLY)
configure_file(mshr-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mshr-config.cmake" @ONLY)

# Install the FooBarConfig.cmake and FooBarConfigVersion.cmake
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mshrConfig.cmake" DESTINATION "${SHARE_DIR}/cmake" COMPONENT dev)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mshr-config.cmake" DESTINATION "${SHARE_DIR}/cmake" COMPONENT dev)
#install(FILES "Use-mshr.cmake" DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)





option(ENABLE_TESTS "Enable testing" OFF)

if(ENABLE_TESTS)
  enable_testing()
  add_subdirectory(test)
else()
  message(STATUS "Testing disabled")
endif()
