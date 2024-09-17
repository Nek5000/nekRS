#
#  ADD_LIBTOOL_LIBRARY -  Mon Jul  8 13:19:42 EDT 2013
#
#  Use this macro like this:
# ADD_LIBTOOL_LIBRARY(project_name 
#   NAME library
#   SRC_LIST  list
#   [NO_LIBTOOL]
#   [BUILD_SHARED_STATIC type]
#   [DEP_LIBS  list]
#   )
#  
#  the first parameter is the project name.
#  NAME is the name of the library to create
#  SRC_LIST is a list of source files to include
#  NO_LIBTOOL disables creation of the .la file
#  BUILD_SHARED_STATIC specifies to build STATIC libraries, SHARED libraries or BOTH
#  DEP_LIBS is a list of dependency libraries to add
#
#
include(CMakeParseArguments)
include(CreateLibtoolFile)

FUNCTION (ADD_LIBTOOL_LIBRARY)
  set(options NO_LIBTOOL)
  set(oneValueArgs NAME BUILD_SHARED_STATIC)
  set(multiValueArgs SRC_LIST DEP_LIBS LINK_LIBS)
  CMAKE_PARSE_ARGUMENTS(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(VERSION_TEXT "")
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version.c")
    file (STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/version.c" VERSION_TEXT REGEX "^static.*")
  endif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version.c")
  if (${VERSION_TEXT} MATCHES ".*${NAME}.*") 
    STRING(REGEX REPLACE ".*Version ([0-9]+.[0-9]+.[0-9]+).*" "\\1" VERSION_TEXT "${VERSION_TEXT}")
    STRING(REGEX REPLACE "([0-9]+)[.][0-9.]*$" "\\1" MAJOR_VERSION_STRING ${VERSION_TEXT})
    STRING(REGEX REPLACE "[0-9]+[.]([0-9]+)[.][0-9.]*$" "\\1" MINOR_VERSION_STRING ${VERSION_TEXT})
    STRING(REGEX REPLACE "[0-9]+[.][0-9]+[.]([0-9]*)$" "\\1" REVISION_STRING ${VERSION_TEXT})
  else (${VERSION_TEXT} MATCHES ".*${NAME}.*") 
    set(VERSION_TEXT "")
  endif (${VERSION_TEXT} MATCHES ".*${NAME}.*") 
  IF( ANDROID)
    set(VERSION_TEXT "")
  ENDIF( ANDROID)
  if ("${ARG_BUILD_SHARED_STATIC}" STREQUAL "") 
    if (DEFINED BUILD_SHARED_STATIC)
      string(TOUPPER ${BUILD_SHARED_STATIC} BUILD_SHARED_STATIC) 
      set (ARG_BUILD_SHARED_STATIC ${BUILD_SHARED_STATIC})
    else (DEFINED BUILD_SHARED_STATIC)
      set (ARG_BUILD_SHARED_STATIC "BOTH")
      set (BUILD_SHARED_STATIC "BOTH" PARENT_SCOPE)
    endif (DEFINED BUILD_SHARED_STATIC)
  endif ("${ARG_BUILD_SHARED_STATIC}" STREQUAL "") 

  IF (${ARG_BUILD_SHARED_STATIC} STREQUAL "BOTH")
    set (BUILD_STATIC TRUE)
    set (BUILD_SHARED TRUE)
  elseif (${ARG_BUILD_SHARED_STATIC} STREQUAL "SHARED")
    set (BUILD_STATIC FALSE)
    set (BUILD_SHARED TRUE)
  elseif (${ARG_BUILD_SHARED_STATIC} STREQUAL "STATIC")
    set (BUILD_SHARED FALSE)
    set (BUILD_STATIC TRUE)
  else ()
    set (ARG_BUILD_SHARED_STATIC "BOTH")
    set (BUILD_STATIC TRUE)
    set (BUILD_SHARED TRUE)
  endif ()


  set (BUILD_LIBTOOL TRUE)
  if (${ARG_NO_LIBTOOL})
    set (BUILD_LIBTOOL FALSE)
  endif (${ARG_NO_LIBTOOL})
  set (STATIC_TARGET_NAME ${ARG_NAME})

  if (${BUILD_SHARED}) 
    add_library( ${ARG_NAME} SHARED ${ARG_SRC_LIST})
    SET_TARGET_PROPERTIES(${ARG_NAME} PROPERTIES LINKER_LANGUAGE C)
    TARGET_LINK_LIBRARIES(${ARG_NAME} ${ARG_LINK_LIBS})
    if (NOT "${VERSION_TEXT}" STREQUAL "")
	set_target_properties(${ARG_NAME}  PROPERTIES VERSION "${VERSION_TEXT}" SOVERSION "${MAJOR_VERSION_STRING}")
	set_target_properties(${ARG_NAME} PROPERTIES LT_VERSION_CURRENT "${MAJOR_VERSION_STRING}")
	set_target_properties(${ARG_NAME} PROPERTIES LT_VERSION_AGE "${MINOR_VERSION_STRING}")
	set_target_properties(${ARG_NAME} PROPERTIES LT_VERSION_REVISION "${REVISION_STRING}")
    endif (NOT "${VERSION_TEXT}" STREQUAL "")
    INSTALL( TARGETS ${ARG_NAME} DESTINATION lib)
    set (STATIC_TARGET_NAME ${ARG_NAME}-static)
  endif (${BUILD_SHARED}) 

  if (${BUILD_STATIC}) 
    add_library(${STATIC_TARGET_NAME} STATIC ${ARG_SRC_LIST})
    # The library target "${ARG_NAME}-static" has a default OUTPUT_NAME of "${ARG_NAME}-static", so change it.
    if (NOT (${STATIC_TARGET_NAME} STREQUAL ${ARG_NAME})) 
      SET_TARGET_PROPERTIES(${ARG_NAME}-static PROPERTIES OUTPUT_NAME "${ARG_NAME}" )
      # Now the library target "foo-static" will be named "foo.lib" with MS tools.
      # This conflicts with the "foo.lib" import library corresponding to "foo.dll",
      # so we add a "lib" prefix (which is default on other platforms anyway):
      SET_TARGET_PROPERTIES(${ARG_NAME}-static PROPERTIES PREFIX "lib" LINKER_LANGUAGE C)
    endif (NOT (${STATIC_TARGET_NAME} STREQUAL ${ARG_NAME})) 
    SET_TARGET_PROPERTIES(${ARG_NAME} PROPERTIES STATIC_LIB "lib${ARG_NAME}.a")
    TARGET_LINK_LIBRARIES(${STATIC_TARGET_NAME} ${ARG_LINK_LIBS})
    INSTALL(TARGETS ${STATIC_TARGET_NAME} DESTINATION lib)
  else (${BUILD_STATIC})
     set (STATIC_TARGET_NAME ${ARG_NAME})
  endif (${BUILD_STATIC}) 

  SET_TARGET_PROPERTIES(${ARG_NAME} PROPERTIES LT_DEPENDENCY_LIBS "${ARG_DEP_LIBS}")
  SET_TARGET_PROPERTIES(${ARG_NAME} PROPERTIES LT_SHOULDNOTLINK "no")
  if (${BUILD_LIBTOOL})
      CREATE_LIBTOOL_FILE(${ARG_NAME} ${STATIC_TARGET_NAME} /lib)
  endif (${BUILD_LIBTOOL})

ENDFUNCTION()
