######################################################
# - Try to find craydrc (http://directory.fsf.org/wiki/Libfabric)
# Once done this will define
#  CrayDRC_FOUND - System has craydrc
#  CrayDRC_INCLUDE_DIRS - The craydrc include directories
#  CrayDRC_LIBRARIES - The libraries needed to use craydrc
#  CrayDRC_DEFINITIONS - The extra CFLAGS needed to use craydrc

######################################################

# This is a bit of a wierd pattern but it allows to bypass pkg-config and
# manually specify library information
if(NOT (PC_CrayDRC_FOUND STREQUAL "IGNORE"))
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    set(_CrayDRC_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
    if(CrayDRC_ROOT)
      list(INSERT CMAKE_PREFIX_PATH 0 "${CrayDRC_ROOT}")
    elseif(NOT ENV{CrayDRC_ROOT} STREQUAL "")
      list(INSERT CMAKE_PREFIX_PATH 0 "$ENV{CrayDRC_ROOT}")
    endif()
    set(PKG_CONFIG_USE_CrayDRC_CMAKE_PREFIX_PATH ON)

    pkg_check_modules(PC_CrayDRC cray-drc)

    set(CMAKE_PREFIX_PATH ${_CrayDRC_CMAKE_PREFIX_PATH})
    unset(_CrayDRC_CMAKE_PREFIX_PATH)

    if(PC_CrayDRC_FOUND)
      if(BUILD_SHARED_LIBS)
        set(_PC_TYPE)
      else()
        set(_PF_TYPE _STATIC)
      endif()
      set(CrayDRC_INCLUDE_DIRS ${PC_CrayDRC${_PC_TYPE}_INCLUDE_DIRS})
      set(CrayDRC_LIBRARIES ${PC_CrayDRC${_PC_TYPE}_LINK_LIBRARIES})
      set(CrayDRC_DEFINITIONS ${PC_CrayDRC${PC_TYPE}_CFLAGS_OTHER})
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CrayDRC DEFAULT_MSG CrayDRC_LIBRARIES)

if(CrayDRC_FOUND)
  if(NOT TARGET craydrc::craydrc)
    add_library(craydrc::craydrc INTERFACE IMPORTED)
    if(CrayDRC_INCLUDE_DIRS)
      set_target_properties(craydrc::craydrc PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CrayDRC_INCLUDE_DIRS}"
      )
    endif()
    if(CrayDRC_DEFINITIONS)
      set_target_properties(craydrc::craydrc PROPERTIES
        INTERFACE_COMPILE_OPTIONS     "${CrayDRC_DEFINITIONS}"
      )
    endif()
    if(CrayDRC_LIBRARIES)
      set_target_properties(craydrc::craydrc PROPERTIES
        INTERFACE_LINK_LIBRARIES      "${CrayDRC_LIBRARIES}"
      )
    endif()
  endif()
endif()
