#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

function(setup_version BASE)
  set(ver_tweak)
  if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/.git)
    if(NOT GIT_COMMAND)
      find_program(GIT_COMMAND git)
    endif()
    if(GIT_COMMAND)
      execute_process(
        COMMAND git
          --git-dir=${CMAKE_CURRENT_LIST_DIR}/.git
          describe --tags --match v${BASE}
        RESULT_VARIABLE res
        OUTPUT_VARIABLE out
        ERROR_QUIET
      )
      if(res EQUAL 0 AND out MATCHES "[^-]*-([^-]*)-g([a-f0-9]*)")
        set(ver_tweak ${CMAKE_MATCH_1})
        set(ver_gitsha ${CMAKE_MATCH_2})
      endif()
    endif()
  endif()

  if(ADIOS2_USE_PIP AND EXISTS "${CMAKE_SOURCE_DIR}/VERSION.TXT")
    file(READ "VERSION.TXT" version_from_file)
    set(ADIOS2_VERSION ${version_from_file} PARENT_SCOPE)
  else()
    if(ver_tweak)
      set(ADIOS2_VERSION ${BASE}.${ver_tweak} PARENT_SCOPE)
    else()
      set(ADIOS2_VERSION ${BASE} PARENT_SCOPE)
    endif()
  endif()

  if(ver_gitsha)
    set(ADIOS2_VERSION_GIT_SHA ${ver_gitsha} PARENT_SCOPE)
  else()
    unset(ADIOS2_VERSION_GIT_SHA PARENT_SCOPE)
  endif()

  set(ADIOS2_LIBRARY_VERSION ${BASE} PARENT_SCOPE)

  string(REGEX MATCH "^([0-9]+\.[0-9]+)" ignore ${BASE})
  set(ADIOS2_LIBRARY_SOVERSION ${CMAKE_MATCH_1} PARENT_SCOPE)
endfunction()

function(adios_option name description default)
  set(ADIOS2_USE_${name} ${default} CACHE STRING "${description}")
  set_property(CACHE ADIOS2_USE_${name} PROPERTY
    STRINGS "ON;TRUE;AUTO;OFF;FALSE"
  )
endfunction()


function(message_pad msg out_len out_msg)
  string(LENGTH "${msg}" msg_len)
  if(NOT (msg_len LESS out_len))
    set(${out_msg} "${msg}" PARENT_SCOPE)
  else()
    math(EXPR pad_len "${out_len} - ${msg_len}")
    string(RANDOM LENGTH ${pad_len} pad)
    string(REGEX REPLACE "." " " pad "${pad}")
    set(${out_msg} "${msg}${pad}" PARENT_SCOPE)
  endif()
endfunction()


function(python_add_test)
  set(options)
  set(oneValueArgs
      NAME
  )
  # EXEC_WRAPPER: Any extra arguments to pass on the command line before test case
  # SCRIPT: Script name and corresponding comand line inputs
  set(multiValueArgs EXEC_WRAPPER SCRIPT)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" "${ARGN}")
  add_test(NAME ${ARGS_NAME}
    COMMAND ${ARGS_EXEC_WRAPPER} $<TARGET_FILE:Python::Interpreter> ${CMAKE_CURRENT_SOURCE_DIR}/${ARGS_SCRIPT}
  )

  if(UNIX)
    set_property(TEST ${ARGS_NAME} PROPERTY
      ENVIRONMENT
        "PYTHONPATH=${ADIOS2_BINARY_DIR}/${CMAKE_INSTALL_PYTHONDIR}:$ENV{PYTHONPATH}"
    )
  else()
    set_property(TEST ${ARGS_NAME} PROPERTY
      ENVIRONMENT
        "PYTHONPATH=${ADIOS2_BINARY_DIR}/${CMAKE_INSTALL_PYTHONDIR};$ENV{PYTHONPATH}"
        "PATH=$<TARGET_FILE_DIR:adios2_py>;$ENV{PATH}"
    )
  endif()
endfunction()


function(GenerateADIOSHeaderConfig)
  set(ADIOS2_CONFIG_DEFINES)
  foreach(OPT IN LISTS ARGN)
    string(TOUPPER ${OPT} OPT_UPPER)
    string(APPEND ADIOS2_CONFIG_DEFINES "
/* CMake Option: ADIOS2_USE_${OPT}=OFF */
#cmakedefine ADIOS2_HAVE_${OPT_UPPER}
")
    if(ADIOS2_HAVE_${OPT})
      set(ADIOS2_HAVE_${OPT_UPPER} 1)
      string(APPEND ADIOS2_CONFIG_FEATURE_LIST "\"${OPT_UPPER}\",")
    else()
      set(ADIOS2_HAVE_${OPT_UPPER})
    endif()
  endforeach()
  string(APPEND ADIOS2_CONFIG_FEATURE_LIST "nullptr")

  configure_file(
    ${ADIOS2_SOURCE_DIR}/source/adios2/common/ADIOSConfig.h.in
    ${ADIOS2_BINARY_DIR}/source/adios2/common/ADIOSConfig.h.in
  )
  configure_file(
    ${ADIOS2_BINARY_DIR}/source/adios2/common/ADIOSConfig.h.in
    ${ADIOS2_BINARY_DIR}/source/adios2/common/ADIOSConfig.h
  )
endfunction()

macro(__adios2_list_cleanup_for_bash var)
  if(${var})
    list(REMOVE_DUPLICATES ${var})
  endif()
  string(REPLACE ";" " " ${var} "${${var}}")
endmacro()

function(__adios2_list_make_link_args var)
  set(prefixes)
  foreach(lib IN LISTS ${var})
    if(lib MATCHES "^/")
      get_filename_component(lib_dir "${lib}" DIRECTORY)
      list(APPEND prefixes "${lib_dir}")
    endif()
  endforeach()

  set(var_new)
  foreach(prefix IN LISTS prefixes)
    list(APPEND var_new "-L${prefix}")
  endforeach()
  foreach(lib IN LISTS ${var})
    if(lib MATCHES "^/.*/?(${CMAKE_SHARED_LIBRARY_PREFIX}|${CMAKE_STATIC_LIBRARY_PREFIX})(.*)(${CMAKE_SHARED_LIBRARY_SUFFIX}|${CMAKE_STATIC_LIBRARY_SUFFIX})")
      list(APPEND var_new "-l${CMAKE_MATCH_2}")
    else()
      list(APPEND var_new "${lib}")
    endif()
  endforeach()

  set(${var} ${var_new} PARENT_SCOPE)
endfunction()


function(adios2_add_thirdparty_target PackageName TargetName)
  add_library(adios2::thirdparty::${PackageName} INTERFACE IMPORTED GLOBAL)
  target_link_libraries(adios2::thirdparty::${PackageName}
    INTERFACE ${TargetName}
  )
endfunction()

# Setup the test dependencies and fixtures for a given pipeline:
function(SetupTestPipeline basename pipeline do_setup)
  # The ideal way to set these up is via test fixtures.  However since those
  # were only available in >= 3.7 then we can get by with DEPENDS.  Since it's
  # all just setting properties anyways though then we go ahead and set both
  # and if >= 3.7 then the fixtures will just get used, otherwise the DEPENDS
  # will take over as a fallback

  if(do_setup)
    add_test(NAME ${basename}.Setup
      COMMAND ${CMAKE_COMMAND}
        -DDIR_NAME=${CMAKE_CURRENT_BINARY_DIR}/${basename}
        -P ${PROJECT_SOURCE_DIR}/cmake/EmptyDir.cmake
    )
    if(pipeline)
      list(INSERT pipeline 0 "Setup")
    else()
      set(pipeline "Setup;")
    endif()
  endif()
  unset(prev_suffix)
  foreach(curr_step IN LISTS pipeline)
    if(curr_step)
      set(curr_suffix .${curr_step})
    else()
      set(curr_suffix "")
    endif()
    if(DEFINED prev_suffix)
      set_property(TEST ${basename}${prev_suffix} APPEND PROPERTY
        FIXTURES_SETUP ${basename}${curr_suffix}
      )
      set_property(TEST ${basename}${curr_suffix} PROPERTY
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${basename}
      )
      set_property(TEST ${basename}${curr_suffix} APPEND PROPERTY
        DEPENDS ${basename}${prev_suffix}
      )
      set_property(TEST ${basename}${curr_suffix} APPEND PROPERTY
        FIXTURES_REQUIRED ${basename}${curr_suffix}
      )
    endif()
    set(prev_suffix "${curr_suffix}")
  endforeach()
endfunction()

macro(adios2_check_fortran_submodules var)
  include(CheckFortranSourceCompiles)
  CHECK_Fortran_SOURCE_COMPILES([[
module foo
  interface bar
    module subroutine bar_integer(x)
      integer, intent(in) :: x
    end subroutine
    module subroutine bar_real(x)
      real, intent(in) :: x
    end subroutine
  end interface
end module
submodule ( foo ) sub
contains
  module subroutine bar_integer(x)
    integer, intent(in) :: x
  end subroutine
  module subroutine bar_real(x)
    real, intent(in) :: x
  end subroutine
end submodule
program main
end program
]] ${var} SRC_EXT F90)
endmacro()

# Set VERSION/SOVERSION of every shared library target in the given directory
# to be the same as the ADIOS VERSION/SOVERSION.  This is important for the
# third-party libraries bundled with ADIOS2.
function(setup_libversion_dir dir)
  get_directory_property(DIR_TARGETS DIRECTORY "${dir}" BUILDSYSTEM_TARGETS)
  foreach(target ${DIR_TARGETS})
    get_target_property(type ${target} TYPE)
    if (${type} STREQUAL "SHARED_LIBRARY")
      set_target_properties(${target} PROPERTIES
          VERSION ${ADIOS2_LIBRARY_VERSION}
          SOVERSION ${ADIOS2_LIBRARY_SOVERSION}
      )
    endif()
  endforeach()
endfunction()
