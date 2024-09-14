message_start_thirdparty()

set(BUILD_GTEST ON CACHE INTERNAL "")
set(BUILD_GMOCK OFF CACHE INTERNAL "")
if(BUILD_SHARED_LIBS)
  set(gtest_force_shared_crt ON CACHE INTERNAL "")
endif()
set(BUILD_SHARED_LIBS OFF)

add_subdirectory(googletest EXCLUDE_FROM_ALL)

message_end_thirdparty()
