#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#
#
# FindDAOS
# -----------
#
# Try to find the DAOS library
#
# This module defines the following variables:
#
#   DAOS_FOUND        - System has DAOS
#   DAOS_INCLUDE_DIRS - The DAOS include directory
#   DAOS_LIBRARIES    - Link these to use DAOS
#
# and the following imported targets:
#   DAOS::DAOS - The core DAOS library
#
# You can also set the following variable to help guide the search:
#   DAOS_ROOT - The install prefix for DAOS containing the
#                     include and lib folders
#                     Note: this can be set as a CMake variable or an
#                           environment variable.  If specified as a CMake
#                           variable, it will override any setting specified
#                           as an environment variable.

if((NOT DAOS_ROOT) AND (NOT (ENV{DAOS_ROOT} STREQUAL "")))
  set(DAOS_ROOT "$ENV{DAOS_ROOT}")
endif()
if(DAOS_ROOT)
  set(DAOS_INCLUDE_OPTS HINTS ${DAOS_ROOT}/include NO_DEFAULT_PATHS)
  set(DAOS_LIBRARY_OPTS
    HINTS ${DAOS_ROOT}/lib ${DAOS_ROOT}/lib64
    NO_DEFAULT_PATHS
    )
endif()

find_path(DAOS_INCLUDE_DIR daos_api.h ${DAOS_INCLUDE_OPTS})
find_library(DAOS_LIBRARY libdaos.so ${DAOS_LIBRARY_OPTS})
find_library(DFS_LIBRARY libdfs.so ${DAOS_LIBRARY_OPTS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DAOS
  FOUND_VAR DAOS_FOUND
  REQUIRED_VARS DAOS_LIBRARY DFS_LIBRARY DAOS_INCLUDE_DIR
)

if(DAOS_FOUND)
  set(DAOS_INCLUDE_DIRS ${DAOS_INCLUDE_DIR})
  set(DAOS_LIBRARIES ${DAOS_LIBRARY} ${DFS_LIBRARY})
  message(STATUS "DAOS Libraries \"${DAOS_LIBRARIES}\"")
  if(NOT TARGET DAOS::DAOS)
    add_library(DAOS::DAOS UNKNOWN IMPORTED)
    set_target_properties(DAOS::DAOS PROPERTIES
      IMPORTED_LOCATION             "${DAOS_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${DAOS_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES      "${DAOS_LIBRARIES}"
    )
  endif()
endif()
