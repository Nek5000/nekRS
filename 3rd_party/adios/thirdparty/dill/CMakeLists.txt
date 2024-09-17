message_start_thirdparty()

set(BUILD_TESTING OFF)
set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_CMAKEDIR}/thirdparty")

set(DILL_QUIET ON CACHE INTERNAL "")
set(DILL_LIBRARY_PREFIX adios2${ADIOS2_LIBRARY_SUFFIX}_ CACHE INTERNAL "")
set(DILL_INSTALL_HEADERS OFF CACHE INTERNAL "")
set(DILL_INSTALL_PKGCONFIG OFF CACHE INTERNAL "")
set(DILL_RUNTIME_COMPONENT adios2_dill-runtime)
set(DILL_LIBRARY_COMPONENT adios2_dill-libraries)
set(DILL_ARCHIVE_COMPONENT adios2_dill-development)
set(DILL_HEADER_COMPONENT adios2_dill-development)

add_subdirectory(dill)
set(dill_DIR ${CMAKE_CURRENT_BINARY_DIR}/dill CACHE INTERNAL "")
setup_libversion_dir(dill)

message_end_thirdparty()
