# Client maintainer: chuck.atkins@kitware.com

message("DEBUG: PATH: $ENV{PATH}")

set(ENV{atl_ROOT} "$ENV{CI_ROOT_DIR}/atl/install")
set(ENV{PATH} "$ENV{CI_ROOT_DIR}/atl/install/bin:$ENV{PATH}")

set(ENV{ffs_ROOT} "$ENV{CI_ROOT_DIR}/ffs/install")
set(ENV{PATH} "$ENV{CI_ROOT_DIR}/ffs/install/bin:$ENV{PATH}")

string(APPEND dashboard_cache "
")

set(CTEST_BUILD_WARNINGS_AS_ERRORS FALSE)

list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")
include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
