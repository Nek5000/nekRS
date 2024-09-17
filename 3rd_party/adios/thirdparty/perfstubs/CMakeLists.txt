message_start_thirdparty()

if (NOT WIN32 AND BUILD_SHARED_LIBS AND ADIOS2_HAVE_Profiling)
  set(PERFSTUBS_USE_TIMERS ON CACHE INTERNAL "")
  add_subdirectory(perfstubs)
else()
  set(PERFSTUBS_USE_TIMERS OFF CACHE INTERNAL "")
endif()

add_library(perfstubs-interface INTERFACE)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/adios2-perfstubs-interface.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/adios2-perfstubs-interface.h"
)

target_include_directories(perfstubs-interface INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

if (PERFSTUBS_USE_TIMERS)
  target_link_libraries(perfstubs-interface INTERFACE perfstubs)
endif()

message_end_thirdparty()
