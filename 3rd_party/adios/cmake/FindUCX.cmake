#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#
#
# FindUCX
# -----------
#
# Try to find the UCX library
#
# This module defines the following variables:
#
#  UCX_FOUND - System has UCX
#  UCX_INCLUDE_DIRS - The UCX include directories
#  UCX_LIBRARIES - The libraries needed to use UCX
#
# and the following imported targets:
#   ucx::ucx - The UCX library target
#
# You can also set the following variable to help guide the search:
#   UCX_ROOT - The install prefix for UCX containing the
#              include and lib folders
#              Note: this can be set as a CMake variable or an
#                    environment variable.  If specified as a CMake
#                    variable, it will override any setting specified
#                    as an environment variable.

# manually specify library information
include(FindPackageHandleStandardArgs)

if(NOT (PC_UCX_FOUND STREQUAL "IGNORE"))
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    set(_UCX_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
    if(UCX_ROOT)
      list(INSERT CMAKE_PREFIX_PATH 0 "${UCX_ROOT}")
    elseif(NOT ENV{UCX_ROOT} STREQUAL "")
      list(INSERT CMAKE_PREFIX_PATH 0 "$ENV{UCX_ROOT}")
    endif()
    set(PKG_CONFIG_USE_UCX_CMAKE_PREFIX_PATH ON)

    pkg_check_modules(PC_UCX ucx)

    set(CMAKE_PREFIX_PATH ${_UCX_CMAKE_PREFIX_PATH})
    unset(_UCX_CMAKE_PREFIX_PATH)

    if(PC_UCX_FOUND)
      if(BUILD_SHARED_LIBS)
        set(_PC_TYPE)
      else()
        set(_PC_TYPE _STATIC)
      endif()
      set(UCX_INCLUDE_DIRS ${PC_UCX${_PC_TYPE}_INCLUDE_DIRS})
      set(UCX_LIBRARIES ${PC_UCX${_PC_TYPE}_LINK_LIBRARIES})
      set(UCX_LIBRARY_DIRS ${PC_UCX${_PC_TYPE}_LIBRARY_DIRS})
      set(UCX_FOUND ${PC_UCX_FOUND})
      set(UCX_VERSION ${PC_UCX_VERSION})
    endif()
  endif()

  # handle the QUIETLY and REQUIRED arguments and set UCX_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(UCX 
      FOUND_VAR UCX_FOUND
      VERSION_VAR UCX_VERSION
      REQUIRED_VARS UCX_LIBRARIES)


else()

  find_path(ucx_ucp_INCLUDE_DIR
     NAMES ucp/api/ucp.h
     PATH_SUFFIXES include
     PATHS ${PC_ucx_INCLUDEDIR} "/opt/ucx/" ${UCX_DIR}
  )
  message (STATUS "ucx_ucp_include dir is ${ucx_ucp_INCLUDE_DIR}")
  find_path(ucx_uct_INCLUDE_DIR
    NAMES uct/api/uct.h
    PATH_SUFFIXES include
    PATHS ${PC_ucx_INCLUDEDIR} "/opt/ucx/include" ${UCX_DIR}
  )
  find_path(ucx_ucs_INCLUDE_DIR
    NAMES ucs/config/global_opts.h
    PATH_SUFFIXES include
    PATHS ${PC_ucx_INCLUDEDIR} "/opt/ucx/include" ${UCX_DIR}
  )

  find_library(ucx_ucp_LIBRARY
    NAMES ucp
    PATH_SUFFIXES lib
    PATHS ${PC_ucx_LIBDIR} "/opt/ucx/lib" ${UCX_DIR}
  )
  find_library(ucx_uct_LIBRARY
    NAMES uct
    PATH_SUFFIXES lib
    PATHS ${PC_ucx_LIBDIR} "/opt/ucx/lib" ${UCX_DIR}
  )
  find_library(ucx_ucs_LIBRARY
    NAMES ucs
    PATH_SUFFIXES lib
    PATHS ${PC_ucx_LIBDIR} "/opt/ucx/lib" ${UCX_DIR}
  )

  set(_ucx_VER_FILE "${ucp_INCLUDE_DIR}/ucp/api/ucp_version.h")
  if(EXISTS "${_ucx_VER_FILE}")
    file(READ "${_ucx_VER_FILE}" _ver)
    string(REGEX MATCH "#define UCP_API_MAJOR *([0-9]*)" _ ${_ver})
    set(_major ${CMAKE_MATCH_1})
    string(REGEX MATCH "#define UCP_API_MINOR *([0-9]*)" _ ${_ver})
    set(_minor ${CMAKE_MATCH_1})
    set(UCX_VERSION "${_major}.${_minor}")
  endif()

  find_package_handle_standard_args(UCX
    FOUND_VAR UCX_FOUND
    REQUIRED_VARS
      ucx_ucp_LIBRARY
      ucx_uct_LIBRARY
      ucx_ucs_LIBRARY
      ucx_ucp_INCLUDE_DIR
      ucx_uct_INCLUDE_DIR
    VERSION_VAR UCX_VERSION
  )

  if(UCX_FOUND)
    set(UCX_LIBRARIES
      ${ucx_ucp_LIBRARY}
      ${ucx_uct_LIBRARY}
      ${ucx_ucs_LIBRARY}
    )
    set(UCX_INCLUDE_DIRS
      ${ucx_ucp_INCLUDE_DIR}
      ${ucx_uct_INCLUDE_DIR}
    )
  endif()
endif()

if(UCX_FOUND)
  if(NOT TARGET ucx::ucx)
    add_library(ucx::ucx INTERFACE IMPORTED)
    if(UCX_INCLUDE_DIRS)
      set_target_properties(ucx::ucx PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UCX_INCLUDE_DIRS}"
      )
    endif()
    if(UCX_LIBRARIES)
      set_target_properties(ucx::ucx PROPERTIES
        INTERFACE_LINK_LIBRARIES      "${UCX_LIBRARIES}"
        INTERFACE_LINK_DIRECTORIES      "${UCX_LIBRARY_DIRS}"
      )
    endif()
  endif()
endif()
