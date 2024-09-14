# Client maintainer: chuck.atkins@kitware.com

set(ENV{CC}  icc)
set(ENV{CXX} icpc)

list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")
include(${CMAKE_CURRENT_LIST_DIR}/unix-common.cmake)
