#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#
#
# PythonModule
# -----------
#
# Try to find the path to a python module
#
# This module defines the following variables:
#
#   PythonModule_${module_NAME}_FOUND - System has the module in question
#
# and the following imported targets:
#   Python::${module_NAME}
#
# This is intented to be called by specifying components for the module name
# and any additional header files and libraries that may be needed to interface
# with it.
#
# For example, say you want to search for mpi4py but you also need it's header
# files to use for interfacing with your own python bindings that use MPI.
# You would call:
#
# find_package(PythonModule REQUIRED COMPONENTS mpi4py mpi4py/mpi4py.h)
#
# If found, this will generate an interface target with the appropriate usage
# requirements.
#
function(__get_filename_label fname var)
  get_filename_component(fname_dir "${fname}" DIRECTORY)
  get_filename_component(fname_nam "${fname}" NAME_WE)
  set(fname_label "${fname_dir}_${fname_nam}")
  string(REGEX REPLACE "[^a-zA-Z0-9_]" "_" fname_label "${fname_label}")
  set(${var} ${fname_label} PARENT_SCOPE)
endfunction()

# Parse out the components
set(module_INCLUDES)
set(module_LIBRARIES)
foreach(comp IN LISTS PythonModule_FIND_COMPONENTS)
  if(comp MATCHES "^[a-zA-Z0-9_]*$")
    set(module_NAME ${comp})
  elseif(comp MATCHES "\\.h$")
    list(APPEND module_INCLUDES "${comp}")
  elseif(comp MATCHES "\\${CMAKE_SHARED_LIBRARY_SUFFIX}")
    list(APPEND module_LIBRARIES "${comp}")
  else()
    message(FATAL_ERROR "${comp}: Unknown python module component")
  endif()
endforeach()
if(NOT module_NAME)
  message(FATAL_ERROR "Unspecified Python module name")
endif()

if(NOT PythonModule_${module_NAME}_FOUND)
  set(PythonModule_${module_NAME}_PATH "" CACHE PATH "Python module ${module_NAME}")
  mark_as_advanced(PythonModule_${module_NAME}_PATH)

  # Try to handle the cross compiling case
  if(NOT PythonModule_${module_NAME}_PATH)
    if(CMAKE_CROSS_COMPILING)
      message(ERROR "find_python_module() invoked in cross-compiling mode, please set the following cache variables appropriately:")
      message("   PythonModule_${module_NAME}_PATH (advanced)")
      message("For details see ${CMAKE_BINARY_DIR}/TryRunResults.cmake")
      file(APPEND "${CMAKE_BINARY_DIR}/TryRunResults.cmake" "
set( PythonModule_${module_NAME}_PATH
  \"PLEASE_FILL_OUT\"
  CACHE PATH \"Python module ${module_NAME}\" FORCE)
")
    else()
      if(NOT Python_Interpreter_FOUND)
        include(CMakeFindDependencyMacro)
        find_dependency(Python COMPONENTS Interpreter)
      endif()
      if(Python_Interpreter_FOUND)
        execute_process(
          COMMAND
            ${Python_EXECUTABLE}
              -c "import ${module_NAME}; print(${module_NAME}.__path__[0])"
          RESULT_VARIABLE result
          OUTPUT_VARIABLE output
          ERROR_VARIABLE error
        )
        if(result EQUAL 0)
          string(STRIP "${output}" output)
          set(PythonModule_${module_NAME}_PATH "${output}"
           CACHE PATH "Python module ${module_NAME}" FORCE)
        endif()
      endif()
    endif()
  endif()

  set(required_vars PythonModule_${module_NAME}_PATH)
  set(include_vars)
  set(library_vars)
  if(PythonModule_${module_NAME}_PATH)
    # Locate interface includes
    foreach(inc IN LISTS module_INCLUDES)
      __get_filename_label("${inc}" inc_label)
      set(inc_var PythonModule_${module_NAME}_${inc_label}_INCLUDE_DIR)
      list(APPEND include_vars ${inc_var})
      find_path(${inc_var} ${inc}
        HINTS ${PythonModule_${module_NAME}_PATH}/include
        NO_CMAKE_FIND_ROOT_PATH NO_DEFAULT_PATH
      )
    endforeach()

    # Locate interface libraries
    foreach(lib IN LISTS module_LIBRARY)
      __get_filename_label("${lib}" lib_label)
      set(lib_var PythonModule_${module_NAME}_${lib_label}_LIBRARY)
      list(APPEND library_vars ${lib_var})
      find_file(${lib_var} ${lib}
        HINTS ${PythonModule_${module_NAME}_PATH}
        NO_CMAKE_FIND_ROOT_PATH NO_DEFAULT_PATH
      )
    endforeach()
    list(APPEND required_vars ${include_vars} ${library_vars})
  endif()

  include(FindPackageHandleStandardArgs)
  set(CMAKE_FIND_PACKAGE_NAME PythonModule_${module_NAME})
  foreach(VAR IN ITEMS REQUIRED QUIETLY VERSION COMPONENTS)
    if(DEFINED PythonModule_FIND_${VAR})
      set(PythonModule_${module_NAME}_FIND_${VAR} "${PythonModule_FIND_${VAR}}")
    else()
      unset(PythonModule_${module_NAME}_FIND_${VAR})
    endif()
  endforeach()

  find_package_handle_standard_args(PythonModule_${module_NAME}
    FOUND_VAR PythonModule_${module_NAME}_FOUND
    REQUIRED_VARS ${required_vars}
  )
endif()

if(PythonModule_${module_NAME}_FOUND AND
  NOT TARGET ${module_NAME})
  add_library(${module_NAME} INTERFACE)
  add_library(Python::${module_NAME} ALIAS ${module_NAME})
  foreach(inc_var IN LISTS include_vars)
    target_include_directories(${module_NAME} SYSTEM INTERFACE ${${inc_var}})
  endforeach()
  foreach(lib_var IN LISTS library_vars)
    target_link_libraries(${module_NAME} INTERFACE ${${lib_var}})
  endforeach()
endif()
