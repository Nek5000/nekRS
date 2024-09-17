cmake_minimum_required(VERSION 3.5)

# The directory label is used for CDash to treat ATL as a subproject of GTKorvo
set(CMAKE_DIRECTORY_LABELS ATL)

project(ATL VERSION 2.2.1 LANGUAGES C)

# Some boilerplate to setup nice output directories
include(GNUInstallDirs)
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/atl
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

if(NOT DEFINED ATL_RUNTIME_COMPONENT)
  set(ATL_RUNTIME_COMPONENT bin)
endif()
if(NOT DEFINED ATL_LIBRARY_COMPONENT)
  set(ATL_LIBRARY_COMPONENT shlib)
endif()
if(NOT DEFINED ATL_ARCHIVE_COMPONENT)
  set(ATL_ARCHIVE_COMPONENT lib)
endif()
if(NOT DEFINED ATL_HEADER_COMPONENT)
  set(ATL_HEADER_COMPONENT dev)
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
      -D_CRT_NONSTDC_NO_DEPRECATE)
  endif()
endif()

# Default to a RelWithDebInfo build if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE RelWithDebInfo)
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

include(CheckTypeSize)
CHECK_TYPE_SIZE(double SIZEOF_DOUBLE)
CHECK_TYPE_SIZE(float SIZEOF_FLOAT)
CHECK_TYPE_SIZE(int SIZEOF_INT)
CHECK_TYPE_SIZE(short SIZEOF_SHORT)

include(CheckIncludeFiles)
CHECK_INCLUDE_FILES(malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H)
CHECK_INCLUDE_FILES(sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES(windows.h HAVE_WINDOWS_H)

include(CheckSymbolExists)
CHECK_SYMBOL_EXISTS(fork unistd.h HAVE_FORK)

if(ATOM_SERVER_HOST)
  set(_atom_server_default "${ATOM_SERVER_HOST}")
else()
  set(_atom_server_default "atomhost.cercs.gatech.edu")
endif()
set(ATL_ATOM_SERVER_HOST "${_atom_server_default}" CACHE STRING
  "Default Atom server")
mark_as_advanced(ATL_ATOM_SERVER_HOST)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

set(ATL_LIBRARY_PREFIX "" CACHE STRING
  "Prefix to prepend to the output library name")
mark_as_advanced(ATL_LIBRARY_PREFIX)

add_library(atl atom.c attr.c lookup3.c tclHash.c)
if(WIN32)
  target_link_libraries(atl wsock32 ws2_32)
endif()
add_library(atl::atl ALIAS atl)
option(ATL_INSTALL_HEADERS "Install ATL header files" ON)
mark_as_advanced(ATL_INSTALL_HEADERS)
set_target_properties(atl PROPERTIES
  OUTPUT_NAME ${ATL_LIBRARY_PREFIX}atl
  VERSION ${ATL_VERSION}
  SOVERSION ${ATL_VERSION_MAJOR})
if(ATL_INSTALL_HEADERS)
  set_target_properties(atl PROPERTIES
    PUBLIC_HEADER atl.h)
endif()
target_include_directories(atl PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Setup pkgconfig
option(ATL_INSTALL_PKGCONFIG "Install ATL pkgconfig files" ON)
mark_as_advanced(ATL_INSTALL_PKGCONFIG)
if(ATL_INSTALL_PKGCONFIG)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/atl.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/atl.pc
    @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/atl.pc
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
    COMPONENT ${ATL_HEADER_COMPONENT})
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/atl-config.in
    ${CMAKE_CURRENT_BINARY_DIR}/atl-config
    @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/atl-config
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    COMPONENT ${ATL_HEADER_COMPONENT})
endif()

if(ATL_INSTALL_HEADERS)
  install(FILES atl.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    COMPONENT ${ATL_HEADER_COMPONENT})
endif()
set(namelink_component_args)
if(NOT CMAKE_VERSION VERSION_LESS 3.12)
  set(namelink_component_args NAMELINK_COMPONENT ${ATL_HEADER_COMPONENT})
endif()
install(TARGETS atl
  # IMPORTANT: Add the foo library to the "export-set"
  EXPORT atl-targets
  RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT ${ATL_RUNTIME_COMPONENT}
  LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT ${ATL_LIBRARY_COMPONENT} ${namelink_component_args}
  ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT ${ATL_ARCHIVE_COMPONENT}
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT ${ATL_HEADER_COMPONENT})


option(ATL_LIBRARIES_ONLY "Don't built the ATL utility executables" OFF)
if(NOT ATL_LIBRARIES_ONLY)
  if(NOT WIN32)
    # Atom server not ported to windows
    add_executable(atom_server atom_server.c)
    target_link_libraries(atom_server atl)
    install(TARGETS atom_server
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT ${ATL_RUNTIME_COMPONENT})
  endif()
  add_executable(attr_dump attr_dump.c)
  target_link_libraries(attr_dump atl)
  if(WIN32)
    target_link_libraries(attr_dump wsock32 ws2_32)
  endif()

  install(TARGETS attr_dump
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT ${ATL_RUNTIME_COMPONENT})
endif()

include(CTest)
mark_as_advanced(BUILD_TESTING)
if(BUILD_TESTING)
  add_executable(attr_test attr_test.c)
  target_link_libraries(attr_test atl)
  add_test(NAME atl_attr_test COMMAND attr_test)

  add_executable(atom_test atom_test.c)
  target_link_libraries(atom_test atl)

  add_executable(atom_check atom_check.c)
  target_link_libraries(atom_check atl)
  add_test(NAME atl_atom_check COMMAND atom_check) 
endif()

# Setup packaging and configs
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/atl-config-version.cmake"
  COMPATIBILITY SameMajorVersion)
configure_file(
  atl-config.cmake.in
  "${PROJECT_BINARY_DIR}/atl-config.cmake"
  @ONLY)
 
# Add all targets to the build-tree export set
export(TARGETS atl NAMESPACE atl::
  FILE "${PROJECT_BINARY_DIR}/atl-targets.cmake")
 
# Install the atl-config.cmake and atl-config-version.cmake
install(
  FILES
    "${PROJECT_BINARY_DIR}/atl-config.cmake"
    "${PROJECT_BINARY_DIR}/atl-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
  COMPONENT ${ATL_HEADER_COMPONENT})
 
# Install the export set for use with the install-tree
install(
  EXPORT atl-targets
  NAMESPACE atl::
  DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
  COMPONENT ${ATL_HEADER_COMPONENT})

set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CMAKE_CURRENT_SOURCE_DIR}/atl.supp)

configure_file(
  ${PROJECT_SOURCE_DIR}/CTestCustom.cmake.in
  ${PROJECT_BINARY_DIR}/CTestCustom.cmake
  @ONLY
)

# Display status message for important variables
option(ATL_QUIET "Suppress info output at the end of configure" OFF)
mark_as_advanced(ATL_QUIET)
if(NOT ATL_QUIET)
  message(STATUS )
  message( STATUS "----------------------------------------------------------------------------" )
  message( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
  message( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
  message( STATUS "ATL_ATOM_SERVER_HOST = ${ATL_ATOM_SERVER_HOST}" )
  message( STATUS "ATL_LIBRARIES_ONLY = ${ATL_LIBRARIES_ONLY}" )
  message( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
  message( STATUS "BUILD_TESTING = ${BUILD_TESTING}" )
  message( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
  message( STATUS "----------------------------------------------------------------------------" )
  message( STATUS )
endif()
