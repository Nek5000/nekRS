# Client maintainer: chuck.atkins@kitware.com

set(dashboard_cache "
")

set(CTEST_CMAKE_GENERATOR "Visual Studio 16")
set(CTEST_CMAKE_GENERATOR_PLATFORM x64)

list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")
include(${CMAKE_CURRENT_LIST_DIR}/windows-common.cmake)
