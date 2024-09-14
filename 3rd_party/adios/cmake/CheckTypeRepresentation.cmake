# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
CheckTypeRepresentation
-----------------------

Check the bit representation of a type.

.. command:: check_float_type_representation

  .. code-block:: cmake

    check_float_type_representation(<type> <var> [LANGUAGE <language>])

  Check if the floating point type exists and determine its bit format.
  ``<type>`` can be any floating point type supported by the compiler, but is
  typically one of ``float``, ``double``, or ``long double``. If the type
  exists, then the variable ``HAVE_<var>`` is set to true, otherwise it is set
  to false. If the type exists, and ``check_float_type_representation()`` was
  able to detect the bit representation of ``<type>``, then the variable
  ``<var>`` is set to its name. If the type exists, but
  ``check_float_type_representation()`` was not able to detect its type, then
  the variable ``<var>`` is set to ``UNKNOWN``. If the type does not exist,
  then ``<var>`` is set to an empty string.

  If ``LANGUAGE`` is set, the specified compiler will be used to perform the
  check. Acceptable values are ``C`` and ``CXX``.

  ``check_float_type_representation()`` understands the following floating
  point formats:

  ``FLOAT_HP_IEEE754_LE``
  ``FLOAT_HP_IEEE754_BE``

    IEEE 754 binary half-precision floating point number. Includes a ``LE`` or
    ``BE`` suffix, depending on whether the number is little endian or big
    endian, respectively.

  ``FLOAT_SP_IEEE754_LE``
  ``FLOAT_SP_IEEE754_BE``

    IEEE 754 binary single-precision floating point number. Includes a ``LE``
    or ``BE`` suffix, depending on whether the number is little endian or big
    endian, respectively.

  ``FLOAT_DP_IEEE754_LE``
  ``FLOAT_DP_IEEE754_BE``

    IEEE 754 binary double-precision floating point number. Includes a ``LE``
    or ``BE`` suffix, depending on whether the number is little endian or big
    endian, respectively.

  ``FLOAT_QP_IEEE754_LE``
  ``FLOAT_QP_IEEE754_BE``

    IEEE 754 binary quadruple-precision floating point number. Includes a
    ``LE`` or ``BE`` suffix, depending on whether the number is little endian
    or big endian, respectively.

  ``FLOAT_QP_IEEE754_80_LE``

    IEEE 754 binary quadruple-precision floating point number with only 80 real
    bits. Typically found on x86 architectures. Since x86 is always a little
    endian architecture, there is no big endian variant of this format.

  ``FLOAT_EP_X86_LE``
  ``FLOAT_EP_X86_64_LE``

    x86 binary 80-bit floating point number. The ``X86`` variant has 16 bits of
    padding (96 bits total), and the ``X86_64`` variant has 48 bits of padding
    (128 bits total). Since x86 is always a little endian architecture, there
    is no big endian variant of this format.

  ``FLOAT_EP_IBM_EXTENDED_LE``
  ``FLOAT_EP_IBM_EXTENDED_BE``

    IBM extended precision floating point number. Typically found on PowerPC
    architectures. Includes a ``LE`` or ``BE`` suffix, depending on whether the
    number is little endian or big endian, respectively.

    This format consists of two IEEE 754 binary double-precision floating point
    numbers back to back. The full floating point number is the sum of these
    two numbers. The first number contains the first half of the significand,
    and the second number contains the second half.
#]=======================================================================]

get_filename_component(__check_type_representation_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)

function(__check_type_representation_impl type value var builtin header_list language)
  # Include header files.
  set(headers)
  if(builtin)
    if(HAVE_SYS_TYPES_H)
      string(APPEND headers "#include <sys/types.h>\n")
    endif()
    if(HAVE_STDINT_H)
      string(APPEND headers "#include <stdint.h>\n")
    endif()
    if(HAVE_STDDEF_H)
      string(APPEND headers "#include <stddef.h>\n")
    endif()
  endif()
  foreach(h ${header_list})
    string(APPEND headers "#include \"${h}\"\n")
  endforeach()

  # Perform the check.

  if(language STREQUAL "C")
    set(cfg_src ${__check_type_representation_dir}/CheckTypeRepresentationCompile.c.in)
    set(src ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CheckTypeRepresentation/${var}_COMPILE.c)
  elseif(language STREQUAL "CXX")
    set(cfg_src ${__check_type_representation_dir}/CheckTypeRepresentationCompile.c.in)
    set(src ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CheckTypeRepresentation/${var}_COMPILE.cpp)
  elseif(language STREQUAL "Fortran")
    set(cfg_src ${__check_type_representation_dir}/CheckTypeRepresentationCompile.f.in)
    set(src ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CheckTypeRepresentation/${var}_COMPILE.f90)
  else()
    message(FATAL_ERROR "Unknown language:\n  ${language}\nSupported languages: C, CXX, Fortran.\n")
  endif()
  set(bin ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CheckTypeRepresentation/${var}_COMPILE.bin)
  configure_file(${cfg_src} ${src} @ONLY)
  try_compile(HAVE_${var} ${CMAKE_BINARY_DIR} ${src}
    COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
    LINK_OPTIONS ${CMAKE_REQUIRED_LINK_OPTIONS}
    LINK_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES}
    CMAKE_FLAGS
      "-DCOMPILE_DEFINITIONS:STRING=${CMAKE_REQUIRED_FLAGS}"
      "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
    OUTPUT_VARIABLE output
    COPY_FILE ${bin}
    )

  if(HAVE_${var})
    # TODO MacOS multi-arch support
    file(READ ${bin} content HEX)

    set(_header "3136206279746520686561646572205b") # string(HEX "16 byte header [" _header)
    set(_footer "5d203332206279746520666f6f746572") # string(HEX "] 32 byte footer" _footer)
    if(content MATCHES "${_header}([0-9a-f]*)${_footer}")
      set(REPR_${var} "${CMAKE_MATCH_1}" CACHE INTERNAL "${type} representation of ${value}")
    else()
      message(SEND_ERROR "Could not find type information while trying to determine representation of ${type}")
    endif()

    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Determining representation of ${type} at compile-time passed with the following output:\n${output}\n\n")
  else()
    # The check failed to compile.
    file(READ ${src} content)
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Determining representation of ${type} at compile-time failed with the following output:\n${output}\n${src}:\n${content}\n\n")

    # Not all compilers support INFINITY and NAN as a const expression.
    # Don't give up hope - try running it instead.
    if(language STREQUAL "C")
      set(cfg_src ${__check_type_representation_dir}/CheckTypeRepresentationRun.c.in)
      set(src ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CheckTypeRepresentation/${var}_RUN.c)
    elseif(language STREQUAL "CXX")
      set(cfg_src ${__check_type_representation_dir}/CheckTypeRepresentationRun.c.in)
      set(src ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CheckTypeRepresentation/${var}_RUN.cpp)
    elseif(language STREQUAL "Fortran")
      set(cfg_src)
    else()
      message(FATAL_ERROR "Unknown language:\n  ${language}\nSupported languages: C, CXX, Fortran.\n")
    endif()
    if(cfg_src)
      configure_file(${cfg_src} ${src} @ONLY)
      try_run(RAN_${var} HAVE_${var} ${CMAKE_BINARY_DIR} ${src}
        COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
        LINK_OPTIONS ${CMAKE_REQUIRED_LINK_OPTIONS}
        LINK_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES}
        CMAKE_FLAGS
          "-DCOMPILE_DEFINITIONS:STRING=${CMAKE_REQUIRED_FLAGS}"
          "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
        COMPILE_OUTPUT_VARIABLE compile_output
        RUN_OUTPUT_VARIABLE run_output
        )
    endif()

    if(NOT HAVE_${var})
      # The check failed to compile.
      if(cfg_src)
        file(READ ${src} content)
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
          "Determining representation of ${type} at run-time failed to compile with the following output:\n${compile_output}\n${src}:\n${content}\n\n")
      endif()
      set(REPR_${var} "" CACHE INTERNAL "${type} representation of ${value}")
    else()
      file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
        "Determining representation of ${type} at run-time compiled with the following output:\n${compile_output}\n\n")
      if(NOT RAN_${var} EQUAL 0)
        # The check failed to run.
        file(READ ${src} content)
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
          "Determining representation of ${type} at run-time failed to run with the following output:\n${run_output}\n${src}:\n${content}\n\n")
        set(REPR_${var} "" CACHE INTERNAL "${type} representation of ${value}")
      else()
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
          "Determining representation of ${type} at run-time ran with the following output:\n${run_output}\n\n")
        set(REPR_${var} "${run_output}" CACHE INTERNAL "${type} representation of ${value}")
      endif()
    endif()
  endif()
endfunction()

function(check_float_type_representation type var)
  if(NOT DEFINED ${var})
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Check representation of ${type}")
    endif()

    cmake_parse_arguments(_cftr "" "LANGUAGE" "" ${ARGN})

    if(NOT DEFINED _cftr_LANGUAGE)
      set(_cftr_LANGUAGE C)
    endif()

    include(Internal/FloatRepresentationTable)

    set(HAVE_${var} TRUE CACHE INTERNAL "check_float_type_representation(${type}): TRUE")

    foreach(_test IN LISTS _tests)
      foreach(_sample IN LISTS _test_${_test})
        __check_type_representation_impl("${type}" "${_sample}" ${var}_${_test} 0 "math.h" ${_cftr_LANGUAGE})
        if(HAVE_${var}_${_test})
          break()
        endif()
      endforeach()
      message(TRACE "HAVE_${var}_${_test}: ${HAVE_${var}_${_test}}")
      message(TRACE "REPR_${var}_${_test}: ${REPR_${var}_${_test}}")

      if(NOT HAVE_${var}_${_test})
        set(HAVE_${var} FALSE CACHE INTERNAL "check_float_type_representation(${type}): FALSE")
      endif()
    endforeach()
    message(TRACE "HAVE_${var}: ${HAVE_${var}}")

    set(_repr_name)
    if(HAVE_${var})
      set(_repr_name "UNKNOWN")
      foreach(_repr IN LISTS _reprs)
        set(_repr_name "${_repr}")
        foreach(_test IN LISTS _tests)
          if(NOT REPR_${var}_${_test} MATCHES "${_repr_${_repr}_${_test}}")
            set(_repr_name "UNKNOWN")
            break()
          endif()
        endforeach()
        if(NOT _repr_name STREQUAL "UNKNOWN")
          break()
        endif()
      endforeach()
    endif()

    if(_repr_name)
      if(NOT CMAKE_REQUIRED_QUIET)
        message(STATUS "Check representation of ${type} - done")
      endif()
      set(${var} "${_repr_name}" CACHE INTERNAL "check_float_type_representation(${type}): ${_repr_name}")
    else()
      if(NOT CMAKE_REQUIRED_QUIET)
        message(STATUS "Check representation of ${type} - failed")
      endif()
      set(${var} "" CACHE INTERNAL "check_float_type_representation(${type}) failed")
    endif()
    message(TRACE "${var}: ${${var}}")
  endif()
endfunction()
