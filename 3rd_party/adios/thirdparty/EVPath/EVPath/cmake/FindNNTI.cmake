######################################################
# Once done this will define
#  NNTI_FOUND - System has nnti
#  NNTI_INCLUDE_DIRS - The nnti include directories
#  NNTI_LIBRARIES - The libraries needed to use nnti

######################################################
set(NNTI_NNTI "" CACHE STRING "Help cmake to find nnti library on your system.")
mark_as_advanced(NNTI_NNTI)
if(NOT NNTI_PREFIX)
  set(NNTI_PREFIX ${NNTI_ROOT})
endif()
if(NOT NNTI_PREFIX)
  set(NNTI_PREFIX $ENV{NNTI_ROOT})
endif()

if(NNTI_PREFIX)
  set(_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
  list(INSERT CMAKE_PREFIX_PATH 0 ${NNTI_PREFIX})
endif()

find_path(NNTI_INCLUDE_DIR Trios_nnti.h)
find_library(NNTI_trios_nnti_LIBRARY trios_nnti)
find_library(NNTI_trios_support_LIBRARY trios_support)
mark_as_advanced(
  NNTI_INCLUDE_DIR NNTI_trios_nnti_LIBRARY NNTI_trios_support_LIBRARY)

if(NNTI_PREFIX)
  set(CMAKE_PREFIX_PATH ${_CMAKE_PREFIX_PATH})
  unset(_CMAKE_PREFIX_PATH)
endif()

######################################################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(nnti DEFAULT_MSG
  NNTI_INCLUDE_DIR NNTI_trios_nnti_LIBRARY NNTI_trios_support_LIBRARY)

if(NNTI_FOUND)
  set(NNTI_LIBRARIES ${NNTI_trios_nnti_LIBRARY} ${NNTI_trios_support_LIBRARY})
  set(NNTI_INCLUDE_DIRS ${NNTI_INCLUDE_DIR})

  if(NOT TARGET nnti::nnti)
    add_library(nnti::nnti UNKNOWN IMPORTED)
    set_target_properties(nnti::nnti PROPERTIES
      IMPORTED_LOCATION "${NNTI_trios_nnti_LIBRARY}"
      INTERFACE_LINK_LIBRARIES "${NNTI_trios_support_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NNTI_INCLUDE_DIR}")
  endif()
endif()
