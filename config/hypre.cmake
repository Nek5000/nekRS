# If user-specified a HYPRE installation, create an imported library target
if (NOT "${HYPRE_DIR}" STREQUAL "")
  if (NOT IS_DIRECTORY ${HYPRE_DIR})
    message(FATAL_ERROR "User specified HYPRE_DIR=${HYPRE_DIR}, but the directory doesn't exist")
  endif()
  message(STATUS "Using HYPRE installation in ${HYPRE_DIR}")
  find_library(HYPRE_FOUND libHYPRE.so PATHS ${HYPRE_DIR}/lib REQUIRED NO_DEFAULT_PATH)
  find_file(HYPRE_H_FOUND HYPRE.h PATHS ${HYPRE_DIR}/include REQUIRED NO_DEFAULT_PATH)
  find_file(HYPRE_PARCSR_LS_H_FOUND HYPRE_parcsr_ls.h PATHS ${HYPRE_DIR}/include REQUIRED NO_DEFAULT_PATH)
  add_library(HYPRE SHARED IMPORTED)
  set_target_properties(HYPRE PROPERTIES
    IMPORTED_LOCATION ${HYPRE_DIR}/lib/libHYPRE.so
    IMPORTED_NO_SONAME FALSE)
  target_include_directories(HYPRE INTERFACE ${HYPRE_DIR}/include)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${HYPRE_DIR}/lib")

# Otherwise, build HYPRE
else()
  set(HYPRE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/hypre)

  # * These two variables are significant to HYPRE's CMakeLists, not our own
  #   HYPRE's CMakeLists leak some variables into parent project, and this is a workaround
  set(HYPRE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE PATH "" FORCE)
  set(HYPRE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "" FORCE)

  set(HYPRE_ENABLE_SINGLE OFF CACHE BOOL "" FORCE)
  set(HYPRE_ENABLE_MIXEDINT ON CACHE BOOL "" FORCE)

  add_subdirectory(${HYPRE_SOURCE_DIR}/src)
  get_property(HYPRE_BINARY_DIR TARGET HYPRE PROPERTY BINARY_DIR)

  # This conflicts with the stdlib "version" header...
  file(REMOVE ${HYPRE_SOURCE_DIR}/src/utilities/version)
endif()
