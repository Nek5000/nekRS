#
#  FIND_CERCS_PROJECT -  Mon Jun 10 12:43:34 EDT 2013
#
#  Use this macro like this:
# FIND_CERCS_PROJECT(project_name 
#   LIBRARY library
#   INCLUDES header1 header2 ...
#   [REQUIRED]
#   [STATIC]
#   [DYNAMIC]
#   [USE_INSTALLED]
#   [VERBOSE]
#   [QUIET]
#   )
#  
#  the first parameter is the project name.
#  LIBRARY is the name of the library to search for (stripped of 'lib'
#     prefix and any postfix)
#  INCLUDES is a list of header files search for
#  REQUIRED fails the build if all not present
#  VERBOSE includes some output about the search path
#  QUIET suppresses the 'found' message
#  USE_INSTALLED avoids searching home and relative-path directories
#
#  the project name is used in directory specs for searching.
#  
#  If both the library and include file are found, then we define the
#  variables:
#  <PROJECT>_FOUND   (where <PROJECT> is the upper-case version of the
#    project argument)
# <PROJECT>_LIB_DIR (suitable for LINK_DIRECTORY calls)
# <PROJECT>_INCLUDE_DIR (suitable for INCLUDE_DIRECTORY calls)
# <PROJECT>_LIBRARIES (full path to a library file)
# HAVE_<include_file>   for each include file found (UPCASED, dot is underscore)
#
include(CMakeParseArguments)

FUNCTION (FIND_CERCS_PROJECT ARG_PROJECT)
  set(options REQUIRED;STATIC;DYNAMIC;USE_INSTALLED;VERBOSE)
  set(oneValueArgs LIBRARY)
  set(multiValueArgs INCLUDES )
  CMAKE_PARSE_ARGUMENTS(ARG "REQUIRED;STATIC;DYNAMIC;USE_INSTALLED;VERBOSE;QUIET" "LIBRARY" "INCLUDES" ${ARGN})

  set (PROJECT_NAME ${ARG_PROJECT})
  set (PROJECT_LIB ${ARG_LIBRARY})
  
  string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPCASE)
  IF (NOT (DEFINED ${PROJECT_NAME_UPCASE}_FOUND))
    if (NOT (DEFINED CercsArch))
      execute_process(COMMAND cercs_arch OUTPUT_VARIABLE CercsArch ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
      MARK_AS_ADVANCED(CercsArch)
    endif()
    set (LIB_SEARCH_PATH)
    set (INC_SEARCH_PATH)
    if (${ARG_STATIC}) 
      set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    endif (${ARG_STATIC}) 
    if (${ARG_DYNAMIC}) 
      set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
    endif (${ARG_DYNAMIC}) 
    if (NOT ${ARG_USE_INSTALLED} )
      if ( NOT ("${CercsArch}" STREQUAL ""))
	list (APPEND LIB_SEARCH_PATH ../${PROJECT_NAME}/${CercsArch} ../../${PROJECT_NAME}/${CercsArch} $ENV{HOME}/${CercsArch}/${PROJECT_NAME}/lib $ENV{HOME}/${CercsArch}/lib )
	list (APPEND INC_SEARCH_PATH ../${PROJECT_NAME}/${CercsArch} ../../${PROJECT_NAME}/${CercsArch} $ENV{HOME}/${CercsArch}/${PROJECT_NAME}/include $ENV{HOME}/${CercsArch}/include )
      endif()
      list (APPEND LIB_SEARCH_PATH ../${PROJECT_NAME}  ../../${PROJECT_NAME} ../${PROJECT_NAME}/build ../../${PROJECT_NAME}/build  $ENV{HOME}/lib  )
      list (APPEND INC_SEARCH_PATH ../${PROJECT_NAME}  ../../${PROJECT_NAME} ../${PROJECT_NAME}/build ../../${PROJECT_NAME}/build  $ENV{HOME}/include  )
    endif (NOT ${ARG_USE_INSTALLED} )
    list (APPEND LIB_SEARCH_PATH ${CMAKE_INSTALL_PREFIX}/lib)
    list (APPEND INC_SEARCH_PATH ${CMAKE_INSTALL_PREFIX}/include)
    IF(EXISTS /users/c/chaos)
      if ( NOT ("${CercsArch}" STREQUAL ""))
	list (APPEND LIB_SEARCH_PATH /users/c/chaos/${CercsArch}/${PROJECT_NAME}/lib /users/c/chaos/${CercsArch}/lib)
	list (APPEND INC_SEARCH_PATH /users/c/chaos/${CercsArch}/${PROJECT_NAME}/include /users/c/chaos/${CercsArch}/include)
      endif()
      list (APPEND LIB_SEARCH_PATH /users/c/chaos/lib )
      list (APPEND INC_SEARCH_PATH /users/c/chaos/include )
    ENDIF(EXISTS /users/c/chaos)
    if (${ARG_VERBOSE})
      message(STATUS "Searching for CERCS PROJECT ${PROJECT_NAME} ${PROJECT_LIB} library in ${LIB_SEARCH_PATH}")
    endif()
    FIND_LIBRARY (${PROJECT_NAME_UPCASE}_LIBRARIES PATHS ${LIB_SEARCH_PATH} NAMES ${PROJECT_LIB} NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
    if ("${${PROJECT_NAME_UPCASE}_LIBRARIES}" STREQUAL "${PROJECT_NAME_UPCASE}_LIBRARIES-NOTFOUND") 
      unset({${PROJECT_NAME_UPCASE}_LIBRARIES)
    endif()
    FIND_LIBRARY (${PROJECT_NAME_UPCASE}_LIBRARIES HINTS NAMES ${PROJECT_LIB})
    if (${ARG_VERBOSE})
      if ("${${PROJECT_NAME_UPCASE}_LIBRARIES}" STREQUAL "${PROJECT_NAME_UPCASE}_LIBRARIES-NOTFOUND") 
	message(STATUS "      -- NOT FOUND")
      else()
	message(STATUS "      -- Found ${${PROJECT_NAME_UPCASE}_LIBRARIES}")
	set (${PROJECT_NAME_UPCASE}_LIBRARIES ${${PROJECT_NAME_UPCASE}_LIBRARIES} PARENT_SCOPE)
      endif()
    endif()
    foreach (INCLUDE ${ARG_INCLUDES})
      STRING(REGEX REPLACE "[./]" "_" INCLUDE_NAME ${INCLUDE})
      STRING(TOUPPER "HAVE_${INCLUDE_NAME}" INCLUDE_VAR)
      if (${ARG_VERBOSE})
	message(STATUS "Searching for ${INCLUDE} in ${SEARCH_PATH}")
      endif()
      FIND_PATH (${INCLUDE_NAME}_INCLUDE_DIR  HINTS ${INC_SEARCH_PATH} NAMES ${INCLUDE} NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
      if ("${${INCLUDE_NAME}_INCLUDE_DIR}" STREQUAL "${INCLUDE_NAME}_INCLUDE_DIR-NOTFOUND") 
	unset({${INCLUDE_NAME}_INCLUDE_DIR)
      endif()
      FIND_PATH (${INCLUDE_NAME}_INCLUDE_DIR  HINTS NAMES ${INCLUDE} )
      if (NOT("${${INCLUDE_NAME}_INCLUDE_DIR}" STREQUAL "${INCLUDE_NAME}_INCLUDE_DIR-NOTFOUND") )
	set (${INCLUDE_VAR} TRUE PARENT_SCOPE)
	set (${PROJECT_NAME_UPCASE}_INCLUDE_DIR ${${INCLUDE_NAME}_INCLUDE_DIR} PARENT_SCOPE)
	set (${PROJECT_NAME_UPCASE}_INCLUDE_DIR ${${INCLUDE_NAME}_INCLUDE_DIR})
	if (${ARG_VERBOSE})
	  message(STATUS "      -- Found in ${${INCLUDE_NAME}_INCLUDE_DIR}")
	endif()
      else()
	if (${ARG_VERBOSE})
	  message(STATUS "      -- NOT FOUND")
	endif()
      endif()
    endforeach(INCLUDE ${ARG_INCLUDES})
    IF (DEFINED ${PROJECT_NAME_UPCASE}_LIBRARIES)
      get_filename_component ( ${PROJECT_NAME_UPCASE}_LIB_DIR ${${PROJECT_NAME_UPCASE}_LIBRARIES} PATH)
      set (${PROJECT_NAME_UPCASE}_LIB_DIR ${${PROJECT_NAME_UPCASE}_LIB_DIR} PARENT_SCOPE)
    ENDIF()
    include(FindPackageHandleStandardArgs)
    set (${PROJECT_NAME}_FIND_REQUIRED ${ARG_REQUIRED})
    set (${PROJECT_NAME}_FIND_QUIET ${ARG_QUIET})
    if ( NOT ("${ARG_INCLUDES}" STREQUAL ""))
      set (INCLUDE_DIR_VAR ${PROJECT_NAME_UPCASE}_INCLUDE_DIR)
    endif()
    find_package_handle_standard_args(${PROJECT_NAME} DEFAULT_MSG
      ${PROJECT_NAME_UPCASE}_LIB_DIR
      ${PROJECT_NAME_UPCASE}_LIBRARIES
      ${INCLUDE_DIR_VAR}
      )
    set(${PROJECT_NAME_UPCASE}_FOUND ${${PROJECT_NAME_UPCASE}_FOUND} PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()
