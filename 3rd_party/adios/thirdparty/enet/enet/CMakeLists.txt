cmake_minimum_required(VERSION 3.5)

project(ENET VERSION 1.3.14 LANGUAGES C)

# Some boilerplate to setup nice output directories
include(GNUInstallDirs)
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/enet
  CACHE STRING "Installation CMake subdirectory")
mark_as_advanced(CMAKE_INSTALL_CMAKEDIR)

if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
endif()
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
endif()

if(NOT DEFINED ENET_RUNTIME_COMPONENT)
  set(ENET_RUNTIME_COMPONENT bin)
endif()
if(NOT DEFINED ENET_LIBRARY_COMPONENT)
  set(ENET_LIBRARY_COMPONENT shlib)
endif()
if(NOT DEFINED ENET_ARCHIVE_COMPONENT)
  set(ENET_ARCHIVE_COMPONENT lib)
endif()
if(NOT DEFINED ENET_HEADER_COMPONENT)
  set(ENET_HEADER_COMPONENT dev)
endif()
if(NOT DEFINED ENET_DOC_COMPONENT)
  set(ENET_DOC_COMPONENT doc)
endif()

# Default to a RelWithDebInfo build if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE RelWithDebInfo)
endif()

if(WIN32)
  # Automagic to do the DLL / LIB song and dance
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

  # Silence MSVC warnings
  if(CMAKE_C_COMPILER_ID MATCHES "MSVC" OR
     CMAKE_C_SIMULATE_ID MATCHES "MSVC")
    add_definitions(
      -D_CRT_SECURE_NO_DEPRECATE
      -D_CRT_SECURE_NO_WARNINGS
      -D_SCL_SECURE_NO_DEPRECATE
      -D_WINSOCK_DEPRECATED_NO_WARNINGS
      -D_CRT_NONSTDC_NO_DEPRECATE
      -DENET_SRC)
    set (MSVC_PERL_FLAGS "-msvc-long")
  endif()
endif()

include(CMakeDependentOption)

# Setup shared library defaults.  If explicitly specified somehow, then default
# to that.  Otherwise base the default on whether or not shared libs are even
# supported.
get_property(SHARED_LIBS_SUPPORTED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
cmake_dependent_option(BUILD_SHARED_LIBS
  "Build shared libraries (so/dylib/dll)." ${SHARED_LIBS_SUPPORTED}
  "SHARED_LIBS_SUPPORTED" OFF
)
mark_as_advanced(BUILD_SHARED_LIBS)

include(CheckFunctionExists)
check_function_exists(getaddrinfo HAS_GETADDRINFO)
check_function_exists(getnameinfo HAS_GETNAMEINFO)
check_function_exists(gethostbyaddr_r HAS_GETHOSTBYADDR_R)
check_function_exists(gethostbyname_r HAS_GETHOSTBYNAME_R)
check_function_exists(poll HAS_POLL)
check_function_exists(fcntl HAS_FCNTL)
check_function_exists(inet_pton HAS_INET_PTON)
check_function_exists(inet_ntop HAS_INET_NTOP)

include(CheckStructHasMember)
check_struct_has_member("struct msghdr" msg_flags sys/socket.h HAS_MSGHDR_FLAGS)

include(CheckCSourceCompiles)
check_c_source_compiles([=[
#include <sys/socket.h>
int main() { return sizeof(socklen_t); } ]=]
  HAS_SOCKLEN_T
)

set(ENET_LIBRARY_PREFIX "" CACHE STRING
  "Prefix to prepend to the output library name")
mark_as_advanced(ENET_LIBRARY_PREFIX)

add_library(enet
  callbacks.c
  compress.c
  host.c
  list.c
  packet.c
  peer.c
  protocol.c
)
set_target_properties(enet PROPERTIES
  OUTPUT_NAME ${ENET_LIBRARY_PREFIX}enet
  VERSION ${ENET_VERSION}
  SOVERSION ${ENET_VERSION_MAJOR}
)
add_library(enet::enet ALIAS enet)

target_compile_definitions(enet PRIVATE
  $<$<BOOL:${HAS_GETADDRINFO}>:HAS_GETADDRINFO>
  $<$<BOOL:${HAS_GETNAMEINFO}>:HAS_GETNAMEINFO>
  $<$<BOOL:${HAS_GETHOSTBYADDR_R}>:HAS_GETHOSTBYADDR_R>
  $<$<BOOL:${HAS_GETHOSTBYNAME_R}>:HAS_GETHOSTBYNAME_R>
  $<$<BOOL:${HAS_POLL}>:HAS_POLL>
  $<$<BOOL:${HAS_FCNTL}>:HAS_FCNTL>
  $<$<BOOL:${HAS_INET_PTON}>:HAS_INET_PTON>
  $<$<BOOL:${HAS_INET_NTOP}>:HAS_INET_NTOP>
  $<$<BOOL:${HAS_MSGHDR_FLAGS}>:HAS_MSGHDR_FLAGS>
  $<$<BOOL:${HAS_SOCKLEN_T}>:HAS_SOCKLEN_T>
)

if(UNIX)
  target_sources(enet PRIVATE unix.c)
elseif(WIN32)
  target_sources(enet PRIVATE win32.c)
  target_link_libraries(enet wsock32 ws2_32 winmm)
endif()
target_include_directories(enet PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Setup pkgconfig
option(ENET_INSTALL_PKGCONFIG "Install enet pkgconfig files" ON)
mark_as_advanced(ENET_INSTALL_PKGCONFIG)
if(ENET_INSTALL_PKGCONFIG)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/enet.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/enet.pc
    @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/enet.pc
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
    COMPONENT ${ENET_HEADER_COMPONENT})
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/enet-config.in
    ${CMAKE_CURRENT_BINARY_DIR}/enet-config
    @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/enet-config
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    COMPONENT ${ENET_HEADER_COMPONENT})
endif()

option(ENET_INSTALL_HEADERS "Install enet header files" ON)
mark_as_advanced(ENET_INSTALL_HEADERS)
if(ENET_INSTALL_HEADERS)
  install(DIRECTORY include/enet
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT ${ENET_HEADER_COMPONENT}
    FILES_MATCHING PATTERN "*.h"
  )
endif()
set(namelink_component_args)
if(NOT CMAKE_VERSION VERSION_LESS 3.12)
  set(namelink_component_args NAMELINK_COMPONENT ${ENET_HEADER_COMPONENT})
endif()
install(TARGETS enet
  EXPORT enet-targets
  RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT ${ENET_RUNTIME_COMPONENT}
  LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT ${ENET_LIBRARY_COMPONENT} ${namelink_component_args}
  ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT ${ENET_ARCHIVE_COMPONENT}
)

# Setup packaging and configs

# Add all targets to the build-tree export set
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/enet-config-version.cmake"
  COMPATIBILITY SameMajorVersion)
configure_file(enet-config.cmake.in
  "${PROJECT_BINARY_DIR}/enet-config.cmake" @ONLY)
export(TARGETS enet NAMESPACE enet::
  FILE "${PROJECT_BINARY_DIR}/enet-targets.cmake")

# Install the enet-config.cmake and enet-config-version.cmake
install(FILES
  "${PROJECT_BINARY_DIR}/enet-config.cmake"
  "${PROJECT_BINARY_DIR}/enet-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
  COMPONENT ${ENET_HEADER_COMPONENT})

# Install the export set for use with the install-tree
install(EXPORT enet-targets NAMESPACE enet::
  DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
  COMPONENT ${ENET_HEADER_COMPONENT})

# Build documentation
find_package(Doxygen)
cmake_dependent_option(BUILD_DOCS
  "Build documentation." OFF 
  "Doxygen_FOUND" OFF)
if(BUILD_DOCS)
  doxygen_add_docs(enet_docs "${CMAKE_CURRENT_SOURCE_DIR}" ALL)
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html" DESTINATION "${CMAKE_INSTALL_DOCDIR}" COMPONENT ${ENET_DOC_COMPONENT})
endif()
