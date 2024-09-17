message_start_thirdparty()

set(BUILD_TESTING OFF)
set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_CMAKEDIR}/thirdparty")

set(ATL_LIBRARIES_ONLY ON CACHE INTERNAL "")
set(ATL_QUIET ON CACHE INTERNAL "")
set(ATL_LIBRARY_PREFIX adios2${ADIOS2_LIBRARY_SUFFIX}_ CACHE INTERNAL "")
set(ATL_INSTALL_HEADERS OFF CACHE INTERNAL "")
set(ATL_INSTALL_PKGCONFIG OFF CACHE INTERNAL "")
set(ATL_RUNTIME_COMPONENT adios2_atl-runtime)
set(ATL_LIBRARY_COMPONENT adios2_atl-libraries)
set(ATL_ARCHIVE_COMPONENT adios2_atl-development)
set(ATL_HEADER_COMPONENT adios2_atl-development)

add_subdirectory(atl)
set(atl_DIR ${CMAKE_CURRENT_BINARY_DIR}/atl CACHE INTERNAL "")
setup_libversion_dir(atl)

message_end_thirdparty()
