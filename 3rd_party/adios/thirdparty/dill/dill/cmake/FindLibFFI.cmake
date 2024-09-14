# - Try to find LibFFI
# Once done this will define
#  LIBFFI_FOUND - System has LibFFI
#  LIBFFI_INCLUDE_DIRS - The LibFFI include directories
#  LIBFFI_LIBRARIES - The libraries needed to use LibFFI
#  LIBFFI_DEFINITIONS - Compiler switches required for using LibFFI

# This is a bit of a wierd pattern but it allows to bypass pkg-config and
# manually specify library information
if(NOT (PC_LIBFFI_FOUND STREQUAL "IGNORE"))
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    set(_DILL_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
    if(LIBFFI_ROOT)
      list(INSERT CMAKE_PREFIX_PATH 0 "${LIBFFI_ROOT}")
    elseif(NOT ENV{LIBFFI_ROOT} STREQUAL "")
      list(INSERT CMAKE_PREFIX_PATH 0 "$ENV{LIBFFI_ROOT}")
    endif()
    set(PKG_CONFIG_USE_DILL_CMAKE_PREFIX_PATH ON)

    pkg_check_modules(PC_LIBFFI libffi)

    set(CMAKE_PREFIX_PATH ${_DILL_CMAKE_PREFIX_PATH})
    unset(_DILL_CMAKE_PREFIX_PATH)

    if(PC_LIBFFI_FOUND)
      if(BUILD_SHARED_LIBS)
        set(_PC_TYPE)
      else()
        set(_PF_TYPE _STATIC)
      endif()
      set(LIBFFI_INCLUDE_DIRS ${PC_LIBFFI${_PC_TYPE}_INCLUDE_DIRS})
      set(LIBFFI_LIBRARIES ${PC_LIBFFI${_PC_TYPE}_LIBRARIES})
      set(LIBFFI_LIBRARY_DIRS ${PC_LIBFFI${_PC_TYPE}_LIBRARY_DIRS})
      set(LIBFFI_DEFINITIONS ${PC_LIBFFI${PC_TYPE}_CFLAGS_OTHER})
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibFFI DEFAULT_MSG LIBFFI_LIBRARIES)

if(LIBFFI_FOUND)
  if(NOT TARGET libffi::libffi)
    add_library(libffi::libffi INTERFACE IMPORTED)
    if(LIBFFI_INCLUDE_DIRS)
      set_target_properties(libffi::libffi PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LIBFFI_INCLUDE_DIRS}"
      )
    endif()
    if(LIBFFI_DEFINITIONS)
      set_target_properties(libffi::libffi PROPERTIES
        INTERFACE_COMPILE_OPTIONS     "${LIBFFI_DEFINITIONS}"
      )
    endif()
    if(LIBFFI_LIBRARIES)
      set_target_properties(libffi::libffi PROPERTIES
        INTERFACE_LINK_LIBRARIES      "${LIBFFI_LIBRARIES}"
        INTERFACE_LINK_DIRECTORIES    "${LIBFFI_LIBRARY_DIRS}"
      )
    endif()
  endif()
endif()
