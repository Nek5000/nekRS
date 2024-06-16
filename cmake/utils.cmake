function (__MPI_find_compiler LANG QUERY_FLAG OUTPUT_VARIABLE)
  separate_arguments(_MPI_COMPILER_WRAPPER_OPTIONS NATIVE_COMMAND "${QUERY_FLAG}")
  set(DUMMYSRC "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx")
  file(WRITE ${DUMMYSRC} "int main() { return 0; }\n")
  execute_process(
    COMMAND ${MPI_${LANG}_COMPILER} ${_MPI_COMPILER_WRAPPER_OPTIONS} ${DUMMYSRC}
    OUTPUT_VARIABLE  WRAPPER_OUTPUT OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE   WRAPPER_ERR ERROR_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE  WRAPPER_RETURN)
  # Some compiler wrappers will yield spurious zero return values, for example
  # Intel MPI tolerates unknown arguments and if the MPI wrappers loads a shared
  # library that has invalid or missing version information there would be warning
  # messages emitted by ld.so in the compiler output. In either case, we'll treat
  # the output as invalid.
  set(WRAPPER_OUTPUT "${WRAPPER_OUTPUT} ${WRAPPER_ERR}")
  if("${WRAPPER_OUTPUT}" MATCHES "undefined reference|unrecognized|need to set|no version information available|command not found")
    set(WRAPPER_RETURN 255)
  endif()
  # Ensure that no error output might be passed upwards.
  if(NOT WRAPPER_RETURN EQUAL 0)
    unset(WRAPPER_OUTPUT)
  else()
    # Strip leading whitespace
    string(REGEX REPLACE "^ +" "" WRAPPER_OUTPUT "${WRAPPER_OUTPUT}")
  endif()

  unset(UNDERLYING_COMPILER)
  if(WRAPPER_OUTPUT)
    separate_arguments(WRAPPER_OUTPUT)
    list(GET WRAPPER_OUTPUT 0 WRAPPER_OUTPUT_0)
    find_program(UNDERLYING_COMPILER ${WRAPPER_OUTPUT_0})
    message("-- Found MPI_UNDERLYING_COMPILER: ${UNDERLYING_COMPILER}")
    if(NOT EXISTS UNDERLYING_COMPILER)
      unset(UNDERLYING_COMPILER)
    endif()
  endif()

  set(${OUTPUT_VARIABLE} "${UNDERLYING_COMPILER}" PARENT_SCOPE)
endfunction()

function (__MPI_underlying_compiler LANG OUTPUT_VARIABLE)
  foreach (flag IN ITEMS "show" "showme" "craype-verbose")
    __MPI_find_compiler("CXX" "-${flag}" COMPILER)
    if(COMPILER) 
      break()
    endif()
  endforeach()

  if(NOT COMPILER)
    message(FATAL_ERROR "Cannot identify underlying compiler used by ${MPI_${LANG}_COMPILER}")
  endif()

  set(${OUTPUT_VARIABLE} "${COMPILER}" PARENT_SCOPE)
endfunction()
