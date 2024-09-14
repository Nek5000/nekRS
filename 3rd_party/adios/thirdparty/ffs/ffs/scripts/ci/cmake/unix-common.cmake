# Client maintainer: chuck.atkins@kitware.com

set(ENV{atl_ROOT}  "$ENV{CI_ROOT_DIR}/atl/install")
set(ENV{dill_ROOT} "$ENV{CI_ROOT_DIR}/dill/install")

string(APPEND dashboard_cache "
FFS_USE_DILL:BOOL=ON
")

if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")
include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
