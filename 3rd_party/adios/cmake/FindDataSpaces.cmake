#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#
#
# FindDATASPACES
# -----------
#
# Try to find the DataSpaces library
#
# This module defines the following variables:
#
#   DATASPACES_FOUND        - System has DataSpaces
#   DATASPACES_INCLUDE_DIRS - The DataSpaces include directory
#   DATASPACES_LIBRARIES    - Link these to use DataSpaces
#
# and the following imported targets:
#   DataSpaces::DataSpaces - The DataSpaces library target
#
# You can also set the following variable to help guide the search:
#   DATASPACES_ROOT - The install prefix for DataSpaces containing the
#              include and lib folders
#              Note: this can be set as a CMake variable or an
#                    environment variable.  If specified as a CMake
#                    variable, it will override any setting specified
#                    as an environment variable.

if(NOT DATASPACES_FOUND)
  if(NOT DATASPACES_ROOT)
    if(NOT ("$ENV{DATASPACES_ROOT}" STREQUAL ""))
      set(DATASPACES_ROOT "$ENV{DATASPACES_ROOT}")
    else()
      find_program(DSPACES_CONF dspaces_config)
      get_filename_component(DATASPACES_ROOT "${DSPACES_CONF}/../.." ABSOLUTE)
    endif()
  endif()
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    set(_DATASPACES_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
    if(DATASPACES_ROOT)
      list(INSERT CMAKE_PREFIX_PATH 0 "${DATASPACES_ROOT}")
    elseif(NOT ENV{DATASPACES_ROOT} STREQUAL "")
      list(INSERT CMAKE_PREFIX_PATH 0 "$ENV{DATASPACES_ROOT}")
    endif()
    set(PKG_CONFIG_USE_DATASPACES_CMAKE_PREFIX_PATH ON)
    pkg_check_modules(PC_DSPACES dspaces)
    set(CMAKE_PREFIX_PATH ${_DATASPACES_CMAKE_PREFIX_PATH})
    unset(_DATASPACES_CMAKE_PREFIX_PATH)
    if(PC_DSPACES_FOUND)
      if(BUILD_SHARED_LIBS)
        set(_PC_TYPE)
      else()
        set(_PC_TYPE _STATIC)
      endif()
      set(DATASPACES_LIBRARIES ${PC_DSPACES${_PC_TYPE}_LINK_LIBRARIES})
      set(DATASPACES_LIBRARY_HINT ${PC_DSPACES${_PC_TYPE}_LIBRARY_DIRS})
      set(DATASPACES_INCLUDE_DIR ${PC_DSPACES${_PC_TYPE}_INCLUDE_DIRS})
      set(DATASPACES_VERSION ${PC_DSPACES_VERSION})
      find_library(DATASPACES_LIBRARY dspaces HINTS ${DATASPACES_LIBRARY_HINT})
      set(HAVE_DSPACES2 TRUE)
    endif()
  endif()
  if(DATASPACES_ROOT AND NOT PC_DSPACES_FOUND)
    find_program(DSPACES_CONF dspaces_config ${DATASPACES_ROOT}/bin)
    if(DSPACES_CONF)
	  execute_process(COMMAND ${DSPACES_CONF} -l
	    RESULT_VARIABLE RESULT_VAR
	    OUTPUT_VARIABLE DSPACES_CONFIG_STRING
	    ERROR_QUIET
	    OUTPUT_STRIP_TRAILING_WHITESPACE)
      string(REPLACE "-L"  " "   LINK_LIBS_ALL   ${DSPACES_CONFIG_STRING})
      string(REPLACE " -l"  " "   LINK_LIBS_B   ${LINK_LIBS_ALL})
      string(REPLACE " "  ";"   LINK_LIBS   ${LINK_LIBS_B})
	  set(DATASPACES_LIBRARIES)
	  set(DATASPACES_LIBRARY_HINT)
	  foreach(LOOP_VAR ${LINK_LIBS})
	    STRING(FIND ${LOOP_VAR} "-u" DEL_FLG)
	    if(("${DEL_FLG}" EQUAL "-1"))
	      STRING(FIND ${LOOP_VAR} "/" HINT_FLG)
		  if(NOT("${HINT_FLG}" EQUAL "-1"))
		    list(APPEND DATASPACES_LIBRARY_HINT ${LOOP_VAR})
		  else()
		    unset(LOCAL_LIBRARY CACHE)
		    unset(LOCAL_LIBRARY-NOTFOUND CACHE)
		    STRING(FIND ${LOOP_VAR} "stdc++" CPP_FLG)
		    if("${CPP_FLG}" EQUAL "-1")
		      find_library(LOCAL_LIBRARY NAMES "${LOOP_VAR}" HINTS ${DATASPACES_LIBRARY_HINT})
		      if(LOCAL_LIBRARY)
			    list(APPEND DATASPACES_LIBRARIES ${LOCAL_LIBRARY})
		      else()
			    list(APPEND DATASPACES_LIBRARIES ${LOOP_VAR})
		      endif()
		    endif()
		  endif()
	    endif()
  	  endforeach()
	  execute_process(COMMAND ${DSPACES_CONF} -v
        RESULT_VARIABLE RESULT_VAR
        OUTPUT_VARIABLE DATASPACES_VERSION
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)	
    endif()
    set(DATASPACES_INCLUDE_OPTS HINTS ${DATASPACES_ROOT}/include)
    find_path(DATASPACES_INCLUDE_DIR dataspaces.h ${DATASPACES_INCLUDE_OPTS})
    find_library(DATASPACES_LIBRARY dspaces HINTS ${DATASPACES_LIBRARY_HINT}) 
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(DataSpaces
    FOUND_VAR DATASPACES_FOUND
    VERSION_VAR DATASPACES_VERSION
    REQUIRED_VARS DATASPACES_VERSION DATASPACES_INCLUDE_DIR 
      DATASPACES_LIBRARIES #DATASPACES_LIBRARY
  )
  if(DATASPACES_FOUND)
    if(DATASPACES_FOUND AND NOT TARGET DataSpaces::DataSpaces)
      add_library(DataSpaces::DataSpaces UNKNOWN IMPORTED)
      set_target_properties(DataSpaces::DataSpaces PROPERTIES
       	IMPORTED_LOCATION             "${DATASPACES_LIBRARY}"
        INTERFACE_LINK_LIBRARIES      "${DATASPACES_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${DATASPACES_INCLUDE_DIR}"
      )
    endif()
  endif()
endif()
