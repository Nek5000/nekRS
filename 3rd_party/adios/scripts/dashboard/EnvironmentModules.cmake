# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#.rst:
# EnvironmentModules
# ------------------
#
# Make environment module commands available to CMake scripts.  This module
# is compatible with both Lua based Lmod and TCL based EnvironmentModules
#
# Module Command
# ^^^^^^^^^^^^^^
#
# This module searched for the module command in the following variable:
#
# ::
#
#   MODULE_COMMAND - The low level module command to use.  Currently supported
#                    are implementations are the Lua based Lmod and TCL based
#                    EnvironmentModules.  The ENV{MODULESHOME} variable,
#                    usually set by the module environment, is used as a hint
#                    to locate the command.
#
# Provided Functions
# ^^^^^^^^^^^^^^^^^^
#
# This module defines the following functions:
#
# ::
#
#   module(...)                 - Execute an arbitry module command
#   module_swap(out_mod in_mod) - Swap out one currently loaded module for
#                                 another
#   module_list(out_var)        - Retrieve the currently loaded modules,
#                                 making the output available as a properly
#                                 formatted CMake ;list variable.
#   module_avail(out_var)       - Retrieve the availabe modules that can be
#                                 loaded, making the output available as a
#                                 properly formatted CMake ;-seperated list
#                                 variable.

# Execute an aribitrary module command.  Usage:
#   module(cmd arg1 ... argN)
#     Process the given command and arguments as if they were passed
#     directly to the module command in your shell environment.
#   module(
#     COMMAND cmd arg1 .. argN
#     [OUTPUT_VARIABLE out_var]
#     [RESULT_VARIABLE ret_var]
#   )
function(module)
  if(NOT MODULE_COMMAND)
    message(ERROR "Failed to process module command.  MODULE_COMMAND not found")
    return()
  endif()

  set(options)
  set(oneValueArgs OUTPUT_VARIABLE RESULT_VARIABLE)
  set(multiValueArgs COMMAND)
  cmake_parse_arguments(MOD_ARGS
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGV}
  )
  if(NOT MOD_ARGS_COMMAND)
    # If no explicit command argument was given, then treat the calling syntax
    # as: module(cmd args...)
    set(exec_cmd ${ARGV})
  else()
    set(exec_cmd ${MOD_ARGS_COMMAND})
  endif()

  if(MOD_ARGS_OUTPUT_VARIABLE)
    set(err_var_args ERROR_VARIABLE err_var)
   endif()

  execute_process(
    COMMAND mktemp -t module.cmake.XXXXXXXXXXXX
    OUTPUT_VARIABLE tempfile_name
  )
  string(STRIP "${tempfile_name}" tempfile_name)

  # If the $MODULESHOME/init/cmake file exists then assume that the CMake
  # "shell" functionality exits
  if(EXISTS "$ENV{MODULESHOME}/init/cmake")
    execute_process(
      COMMAND ${MODULE_COMMAND} cmake ${exec_cmd}
      OUTPUT_FILE ${tempfile_name}
      ${err_var_args}
      RESULT_VARIABLE ret_var
    )
  else() # fallback to the sh shell and manually convert to CMake
    execute_process(
      COMMAND ${MODULE_COMMAND} sh ${exec_cmd}
      OUTPUT_VARIABLE out_var
      ${err_var_args}
      RESULT_VARIABLE ret_var
    )
  endif()

  # If we executed successfully then process and cleanup the temp file
  if("${ret_var}" EQUAL 0)
    # No CMake shell so we need to process the sh output into CMake code
    if(NOT EXISTS "$ENV{MODULESHOME}/init/cmake")
      file(WRITE ${tempfile_name} "")
      string(REPLACE "\n" ";" out_var "${out_var}")
      foreach(sh_cmd IN LISTS out_var)
        if(sh_cmd MATCHES "^ *unset *([^ ]*)")
          set(cmake_cmd "unset(ENV{${CMAKE_MATCH_1}})")
        elseif(sh_cmd MATCHES "^ *export *([^ ]*)")
          set(cmake_cmd "set(ENV{${CMAKE_MATCH_1}} \"\${${CMAKE_MATCH_1}}\")")
        elseif(sh_cmd MATCHES " *([^ =]*) *= *(.*)")
          set(var_name "${CMAKE_MATCH_1}")
          set(var_value "${CMAKE_MATCH_2}")
          if(var_value MATCHES "^\"(.*[^\\])\"")
            # If it's in quotes, take the value as is
            set(var_value "${CMAKE_MATCH_1}")
          else()
            # Otherwise, strip trailing spaces
            string(REGEX REPLACE "([^\\])? +$" "\\1" var_value "${var_value}")
          endif()
          string(REPLACE "\\ " " " var_value "${var_value}")
          set(cmake_cmd "set(${var_name} \"${var_value}\")")
        else()
          continue()
        endif()
        file(APPEND ${tempfile_name} "${cmake_cmd}\n")
      endforeach()
    endif()

    # Process the change in environment variables
    include(${tempfile_name})
    file(REMOVE ${tempfile_name})
  endif()

  # Push the output back out to the calling scope
  if(MOD_ARGS_OUTPUT_VARIABLE)
    set(${MOD_ARGS_OUTPUT_VARIABLE} "${err_var}" PARENT_SCOPE)
  endif()
  if(MOD_ARGS_RESULT_VARIABLE)
    set(${MOD_ARGS_RESULT_VARIABLE} ${ret_var} PARENT_SCOPE)
  endif()
endfunction(module)

# Swap one module for another
function(module_swap out_mod in_mod)
  module(COMMAND -t swap ${out_mod} ${in_mod} OUTPUT_VARIABLE tmp_out)
endfunction()

# Retrieve the currently loaded modules
function(module_list out_var)
  cmake_policy(SET CMP0007 NEW)
  module(COMMAND -t list OUTPUT_VARIABLE tmp_out)

  # Convert output into a CMake list
  string(REPLACE "\n" ";" ${out_var} "${tmp_out}")

  # Remove title headers and empty entries
  list(REMOVE_ITEM ${out_var} "No modules loaded")
  if(${out_var})
    list(FILTER ${out_var} EXCLUDE REGEX "^(.*:)?$")
  endif()
  list(FILTER ${out_var} EXCLUDE REGEX "^(.*:)?$")

  set(${out_var} ${${out_var}} PARENT_SCOPE)
endfunction()

# Retrieve the list of available modules
function(module_avail out_var)
  cmake_policy(SET CMP0007 NEW)
  module(COMMAND -t avail OUTPUT_VARIABLE tmp_out)

  # Convert output into a CMake list
  string(REPLACE "\n" ";" tmp_out "${tmp_out}")

  set(${out_var})
  foreach(MOD IN LISTS tmp_out)
    # Remove directory entries and empty values
    if(MOD MATCHES "^(.*:)?$")
      continue()
    endif()

    # Convert default modules
    if(MOD MATCHES "^(.*)/$" ) # "foo/"
      list(APPEND ${out_var} ${CMAKE_MATCH_1})
    elseif(MOD MATCHES "^((.*)/.*)\\(default\\)$") # "foo/1.2.3(default)"
      list(APPEND ${out_var} ${CMAKE_MATCH_2})
      list(APPEND ${out_var} ${CMAKE_MATCH_1})
    else()
      list(APPEND ${out_var} ${MOD})
    endif()
  endforeach()

  set(${out_var} ${${out_var}} PARENT_SCOPE)
endfunction()

# Make sure our CMake is new enough
if(CMAKE_VERSION VERSION_LESS 3.6)
  message(FATAL_ERROR
    "The EnvironmentModules interface requires at least CMake v3.6"
  )
endif()

# Make sure we know where the underlying module command is
find_program(MODULE_COMMAND
  NAMES lmod modulecmd
  HINTS ENV MODULESHOME
  PATH_SUFFIXES libexec
)
