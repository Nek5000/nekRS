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

