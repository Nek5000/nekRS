message_start_thirdparty()

set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_CMAKEDIR}/thirdparty")
set(ENET_LIBRARY_PREFIX adios2${ADIOS2_LIBRARY_SUFFIX}_ CACHE INTERNAL "")
set(ENET_INSTALL_HEADERS OFF CACHE INTERNAL "")
set(ENET_INSTALL_PKGCONFIG OFF CACHE INTERNAL "")
set(ENET_RUNTIME_COMPONENT adios2_enet-runtime)
set(ENET_LIBRARY_COMPONENT adios2_enet-libraries)
set(ENET_ARCHIVE_COMPONENT adios2_enet-development)
set(ENET_HEADER_COMPONENT adios2_enet-development)
set(ENET_INSTALL_PKGCONFIG OFF)
set(BUILD_DOCS OFF CACHE BOOL "")

add_subdirectory(enet)
set(enet_DIR ${CMAKE_CURRENT_BINARY_DIR}/enet CACHE INTERNAL "")
setup_libversion_dir(enet)

message_end_thirdparty()
