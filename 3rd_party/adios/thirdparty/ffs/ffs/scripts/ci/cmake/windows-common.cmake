# Client maintainer: chuck.atkins@kitware.com

string(APPEND dashboard_cache "
")

list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")

# the two lines below shouldn't be necessary, but something wrong with vcpkg maybe
set(ENV{dill_DIR} "${CMAKE_CURRENT_LIST_DIR}/../../../../dill/install/lib/cmake/dill")
set(ENV{atl_DIR} "${CMAKE_CURRENT_LIST_DIR}/../../../../atl/install/lib/cmake/atl")

include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
