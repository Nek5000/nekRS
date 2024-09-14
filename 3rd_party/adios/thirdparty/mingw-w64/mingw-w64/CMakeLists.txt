add_library(getopt-obj OBJECT mingw-w64-crt/misc/getopt.c)
target_include_directories(getopt-obj PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/mingw-w64-headers/crt
)

add_library(getopt INTERFACE)
target_include_directories(getopt INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/mingw-w64-headers/crt
)
target_sources(getopt INTERFACE $<TARGET_OBJECTS:getopt-obj>)
add_library(mingw-w64::getopt ALIAS getopt)
