# Client maintainer: chuck.atkins@kitware.com

set(ENV{CC}  gcc)
set(ENV{CXX} g++)
set(ENV{FC}  gfortran)
set(UBSAN_FLAGS "-fsanitize=undefined -fno-sanitize-recover=all -pthread")
set(ENV{CFLAGS}   "${UBSAN_FLAGS}")
set(ENV{CXXFLAGS} "${UBSAN_FLAGS}")
set(ENV{FFLAGS}   "${UBSAN_FLAGS}")
set(ENV{UBSAN_OPTIONS} "print_stacktrace=1")

find_package(EnvModules REQUIRED)
env_module(load mpi)

set(dashboard_cache "
BUILD_TESTING:BOOL=ON
ADIOS2_BUILD_EXAMPLES:BOOL=ON

ADIOS2_USE_HDF5:STRING=ON
ADIOS2_USE_MPI:STRING=ON

HDF5_C_COMPILER_EXECUTABLE:FILEPATH=/usr/lib64/openmpi/bin/h5pcc
HDF5_DIFF_EXECUTABLE:FILEPATH=/usr/lib64/openmpi/bin/h5diff   

MPIEXEC_EXTRA_FLAGS:STRING=--allow-run-as-root --oversubscribe
MPIEXEC_MAX_NUMPROCS:STRING=4
")

set(dashboard_track "Analysis")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_FLAGS "-k -j2")
set(CTEST_MEMORYCHECK_TYPE "UndefinedBehaviorSanitizer")

list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")
include(${CMAKE_CURRENT_LIST_DIR}/ci-common.cmake)
