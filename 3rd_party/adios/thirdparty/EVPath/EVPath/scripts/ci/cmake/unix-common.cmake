# Client maintainer: chuck.atkins@kitware.com

set(ENV{atl_ROOT}  "$ENV{CI_ROOT_DIR}/atl/install")
set(ENV{LD_LIBRARY_PATH} "$ENV{atl_ROOT}/lib:$ENV{LD_LIBRARY_PATH}")

set(ENV{dill_ROOT} "$ENV{CI_ROOT_DIR}/dill/install")
set(ENV{LD_LIBRARY_PATH} "$ENV{dill_ROOT}/lib:$ENV{LD_LIBRARY_PATH}")

set(ENV{ffs_ROOT}  "$ENV{CI_ROOT_DIR}/ffs/install")
set(ENV{LD_LIBRARY_PATH} "$ENV{ffs_ROOT}/lib:$ENV{LD_LIBRARY_PATH}")


string(APPEND dashboard_cache "
")

if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

set(CTEST_BUILD_WARNINGS_AS_ERRORS FALSE)

list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")
include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
