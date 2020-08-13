set(HYPRE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/hypre)

# * These two variables are significant to HYPRE's CMakeLists, not our own
#   HYPRE's CMakeLists leak some variables into parent project, and this is a workaround
set(HYPRE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE PATH "" FORCE)
set(HYPRE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "" FORCE)

# Defaults to OFF in HYPRE's CMakeLists.  This overrides the default.
set(HYPRE_ENABLE_SINGLE ON CACHE BOOL "Use float for HYPRE_Real")

if (EXISTS ${HYPRE_SOURCE_DIR})
  message(STATUS "Using HYPRE source in ${HYPRE_SOURCE_DIR}")
  add_subdirectory(${HYPRE_SOURCE_DIR}/src)
  get_property(HYPRE_BINARY_DIR TARGET HYPRE PROPERTY BINARY_DIR)
else()

  FetchContent_Declare(
    hypre_content
    URL https://github.com/hypre-space/hypre/archive/v${HYPRE_VER}.tar.gz )
  FetchContent_GetProperties(hypre_content)
  if (NOT hypre_content_POPULATED)
    FetchContent_Populate(hypre_content)
  endif()
  set(HYPRE_SOURCE_DIR ${hypre_content_SOURCE_DIR})
  set(HYPRE_BINARY_DIR ${hypre_content_BINARY_DIR})

  # * Exclude from all since HYPRE CMakeLists adds a bunch of targets we don't need
  #   libHYPRE will be build just fine, since we've explicitly declared it as a dependency
  add_subdirectory(${HYPRE_SOURCE_DIR}/src ${HYPRE_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# This conflicts with the stdlib "version" header...
file(REMOVE ${HYPRE_SOURCE_DIR}/src/utilities/version)

