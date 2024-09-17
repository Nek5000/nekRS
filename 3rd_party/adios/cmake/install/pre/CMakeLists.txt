install(CODE "
  message(\"Pre-installation cleanup of CMake files\")
  file(REMOVE_RECURSE \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_CMAKEDIR}\")
"
  COMPONENT adios2_core-development
)
file(GLOB ADIOS2_MODULE_FILES
  "${ADIOS2_SOURCE_DIR}/cmake/Find*.cmake"
  "${ADIOS2_SOURCE_DIR}/cmake/CMake*.cmake"
)
install(FILES ${ADIOS2_MODULE_FILES}
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT adios2_core-development
)
install(DIRECTORY ${ADIOS2_SOURCE_DIR}/cmake/upstream/
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}/upstream
  COMPONENT adios2_core-development
)

