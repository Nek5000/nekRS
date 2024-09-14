# Client maintainer: vicente.bolea@kitware.com
set(ENV{CC}  clang)
set(ENV{CXX} clang++)
set(ENV{FC}  gfortran-11)

set(dashboard_cache "
BUILD_TESTING:BOOL=ON
ADIOS2_BUILD_EXAMPLES:BOOL=ON

ADIOS2_USE_Fortran:BOOL=ON
ADIOS2_USE_MPI:BOOL=OFF
ADISO2_USE_Python:BOOL=ON

CMAKE_C_COMPILER_LAUNCHER=ccache
CMAKE_CXX_COMPILER_LAUNCHER=ccache
CMAKE_C_FLAGS:STRING=-Wall
CMAKE_CXX_FLAGS:STRING=-Wall
CMAKE_Fortran_FLAGS:STRING=-Wall
")

set(ENV{MACOSX_DEPLOYMENT_TARGET} "11.3")
set(CTEST_CMAKE_GENERATOR "Ninja")
list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")
set(CTEST_TEST_ARGS
  # Unclear why these tests currently die.  Disabling until it can be addressed.
  EXCLUDE "Install.Make.*"
)
include(${CMAKE_CURRENT_LIST_DIR}/ci-common.cmake)
