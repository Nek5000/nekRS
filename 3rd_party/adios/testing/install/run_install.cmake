#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

unset(ENV{DESTDIR})
file(REMOVE_RECURSE "${ADIOS2_BINARY_DIR}/testing/install/install")
execute_process(COMMAND "${CMAKE_COMMAND}"
  "-DCMAKE_INSTALL_PREFIX=${ADIOS2_BINARY_DIR}/testing/install/install"
  "-DBUILD_TYPE=${BUILD_TYPE}"
  -P "${ADIOS2_BINARY_DIR}/cmake_install.cmake"
  RESULT_VARIABLE result
  )
if(result)
  message(FATAL_ERROR "Result of installation was ${result}, should be 0")
endif()
