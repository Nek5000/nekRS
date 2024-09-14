configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/perfstubs_api/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/perfstubs_api/config.h"
)

add_library(perfstubs perfstubs_api/timer.c)

set_target_properties(perfstubs PROPERTIES
  OUTPUT_NAME adios2${ADIOS2_LIBRARY_SUFFIX}_perfstubs
  VERSION ${ADIOS2_LIBRARY_VERSION}
  SOVERSION ${ADIOS2_LIBRARY_SOVERSION}
)

target_link_libraries(perfstubs INTERFACE dl m)

target_include_directories(perfstubs PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
install(TARGETS perfstubs EXPORT adios2Exports)
