#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

# Cleanup any lefter old CMake packaging files
file(GLOB LEFTOVER_CONFIGS "${ADIOS2_BINARY_DIR}/ADIOS2*.cmake")
if(LEFTOVER_CONFIGS)
  file(REMOVE_RECURSE ${LEFTOVER_CONFIGS})
endif()

include(CMakePackageConfigHelpers)

# Build interface configs
write_basic_package_version_file(
  ${ADIOS2_BINARY_DIR}/adios2-config-version.cmake
  COMPATIBILITY AnyNewerVersion
)
export(EXPORT adios2Exports
  FILE ${ADIOS2_BINARY_DIR}/adios2-targets.cmake
  NAMESPACE adios2::
)
configure_file(
  ${ADIOS2_SOURCE_DIR}/cmake/adios2-config-common.cmake.in
  ${ADIOS2_BINARY_DIR}/adios2-config-common.cmake
  @ONLY
)
configure_file(
  ${ADIOS2_SOURCE_DIR}/cmake/adios2-config.cmake.in
  ${ADIOS2_BINARY_DIR}/adios2-config.cmake
  @ONLY
)

# Install interface configs
install(
  FILES
    ${ADIOS2_BINARY_DIR}/adios2-config-version.cmake
    ${ADIOS2_BINARY_DIR}/adios2-config-common.cmake
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT adios2_core-development
)
install(EXPORT adios2Exports
  FILE adios2-targets.cmake
  NAMESPACE adios2::
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT adios2_core-development
)
configure_file(
  ${ADIOS2_SOURCE_DIR}/cmake/adios2-config-install.cmake.in
  ${ADIOS2_BINARY_DIR}/adios2-config-install.cmake
  @ONLY
)
install(FILES ${ADIOS2_BINARY_DIR}/adios2-config-install.cmake
  RENAME adios2-config.cmake
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT adios2_core-development
)
