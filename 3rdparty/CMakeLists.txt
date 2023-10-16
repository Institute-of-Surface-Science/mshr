if (NOT USE_SYSTEM_CGAL)
  # Enable ExternalProject CMake module
  include(ExternalProject)
  
  # Set default ExternalProject root directory
  set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/3rdparty)

  message(STATUS "Building CGAL")
  
  message(STATUS "Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
  message(STATUS "Binary dir: ${CMAKE_BINARY_DIR}")
  
  set(CGAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/CGAL-installdir)
  
  if ( CMAKE_BUILD_TYPE )
    set(CGAL_BUILD_TYPE ${CMAKE_BUILD_TYPE})
  else()
    set(CGAL_BUILD_TYPE Debug)
  endif()
  
  # CGAL doesn't allow Developer, RelWithDebInfo and MinSizeRel as build type.
  if ( "${CMAKE_BUILD_TYPE}" STREQUAL "Developer" )
    set(CGAL_BUILD_TYPE "Debug")
  endif()
  if ( "${CMAKE_BUILD_TYPE}" STREQUAL "MinSizeRel"
      OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo" )
    set(CGAL_BUILD_TYPE "Release")
  endif()
  
  message(STATUS "CGAL build type: ${CGAL_BUILD_TYPE}")
  
  # Note: Need to add fPIC in order to be able to link staticly with CGAL
  # FIXME: The CMAKE_CXX_FLAGS=-std=c++11 is a hack to remove CGAL's dependency on boost.thread (see setup_Boost.cmake)
  # This should be removed when CGAL's boost dependency is eventually removed.
  ExternalProject_Add(CGAL
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CGAL
    INSTALL_DIR ${CGAL_INSTALL_DIR}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CGAL_INSTALL_DIR} -DCMAKE_BUILD_TYPE=${CGAL_BUILD_TYPE} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_CXX_FLAGS=-std=c++11 -DCGAL_CXX_FLAGS=-fPIC -DWITH_CGAL_Core:BOOL=OFF -DWITH_CGAL_ImageIO:BOOL=OFF -DWITH_CGAL_Qt5:BOOL=OFF -DWITH_demos:BOOL=OFF -DWITH_examples:BOOL=OFF -DWITH_LEDA:BOOL=OFF -DWITH_GMP:BOOL=ON -DWITH_MPFR:BOOL=ON -DWITH_Eigen3:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=OFF -DBOOST_ROOT=${BOOST_ROOT} -DGMP_INCLUDE_DIR=${GMP_INCLUDE_DIR} -DGMP_LIBRARIES=${GMP_LIBRARIES} -DMPFR_INCLUDE_DIR=${MPFR_INCLUDE_DIR} -DMPFR_LIBRARIES=${MPFR_LIBRARIES} -DEIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR} -DCGAL_HEADER_ONLY:BOOL=${CGAL_HEADER_ONLY}
    
    # -DWITH_ESBTL:BOOL=OFF        \
    # -DWITH_GMPXX:BOOL=OFF        \
    # -DWITH_IPE:BOOL=OFF          \
    # -DWITH_LAPACK:BOOL=OFF       \
    # -DWITH_LEDA:BOOL=OFF         \
    # -DWITH_MPFI:BOOL=OFF         \
    # -DWITH_NTL:BOOL=OFF          \
    # -DWITH_OpenGL:BOOL=OFF       \
    # -DWITH_OpenNL:BOOL=OFF       \
    # -DWITH_QGLViewer:BOOL=OFF    \
    # -DWITH_RS:BOOL=OFF           \
    # -DWITH_RS3:BOOL=OFF          \
    # -DWITH_TAUCS:BOOL=OFF        \
    # -DWITH_ZLIB:BOOL=OFF         \
    
    UPDATE_COMMAND ""
    )
  
  set(CGAL_INCLUDE_DIR ${CGAL_INSTALL_DIR}/include PARENT_SCOPE)
  set(CGAL_LIB_DIR ${CGAL_INSTALL_DIR}/lib PARENT_SCOPE)
  
  # LINK_DIRECTORIES( ${CGAL_INSTALL_DIR}/lib )
  
  # ADD_LIBRARY(CGAL STATIC IMPORTED)
  # SET_TARGET_PROPERTIES(TARGET cgal PROPERTY
  #     IMPORTED_LOCATION ${CGAL_INSTALL_DIR}/lib/libCGAL.a)
  
  # ADD_LIBRARY(CGAL_Core STATIC IMPORTED)
  # SET_TARGET_PROPERTIES(TARGET cgal PROPERTY
  #     IMPORTED_LOCATION ${CGAL_INSTALL_DIR}/lib/libCGAL_Core.a)
endif()

if (NOT USE_SYSTEM_TETGEN)
  add_subdirectory( tetgen1.5.0 )
endif()

set(EXTERNAL_INCLUDE_DIRS "${CGAL_INSTALL_DIR}/include;${TETGEN_INCLUDE_DIR}" PARENT_SCOPE)
set(EXTERNAL_LIBS "tet" PARENT_SCOPE)
set(EXTERNAL_DEFINITIONS "${TETGEN_DEFINITIONS}" PARENT_SCOPE)
