cmake_minimum_required(VERSION 3.5)

# The directory label is used for CDash to treat FFS as a subproject of GTKorvo
set(CMAKE_DIRECTORY_LABELS FFS)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project(FFS VERSION 3.0.0)

# Some boilerplate to setup nice output directories
include(GNUInstallDirs)
if(WIN32 AND NOT CYGWIN)
  set(CMAKE_INSTALL_CMAKEDIR CMake
    CACHE STRING "Installation CMake subdirectory")
else()
  set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ffs
    CACHE STRING "Installation CMake subdirectory")
endif()
mark_as_advanced(CMAKE_INSTALL_CMAKEDIR)

list(INSERT CMAKE_PREFIX_PATH 0 ${CMAKE_INSTALL_PREFIX})
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)
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

list(INSERT CMAKE_PREFIX_PATH 0 ${CMAKE_INSTALL_PREFIX})

if(NOT DEFINED FFS_RUNTIME_COMPONENT)
  set(FFS_RUNTIME_COMPONENT bin)
endif()
if(NOT DEFINED FFS_LIBRARY_COMPONENT)
  set(FFS_LIBRARY_COMPONENT shlib)
endif()
if(NOT DEFINED FFS_ARCHIVE_COMPONENT)
  set(FFS_ARCHIVE_COMPONENT lib)
endif()
if(NOT DEFINED FFS_HEADER_COMPONENT)
  set(FFS_HEADER_COMPONENT dev)
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
      -D_CRT_NONSTDC_NO_DEPRECATE)
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
  "SHARED_LIBS_SUPPORTED" OFF)
mark_as_advanced(BUILD_SHARED_LIBS)
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS -DFFS_SRC)
endif()

# Default to a RelWithDebInfo build if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE RelWithDebInfo)
endif()

set(TARGET_DEP_INC)
set(TARGET_DEP_LIBS)
set(PKG_DEP_PKG)
set(PKG_DEP_LIBS)

include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckTypeSize)
include(CheckStructHasMember)
include(CMakeFindDependencyMacro)
include(TestBigEndian)

CHECK_TYPE_SIZE("double" SIZEOF_DOUBLE)
CHECK_TYPE_SIZE("float" SIZEOF_FLOAT)
CHECK_TYPE_SIZE("int" SIZEOF_INT)
CHECK_TYPE_SIZE("off_t" SIZEOF_OFF_T)
CHECK_TYPE_SIZE("long" SIZEOF_LONG)
CHECK_TYPE_SIZE("long double" SIZEOF_LONG_DOUBLE)
CHECK_TYPE_SIZE("long long" SIZEOF_LONG_LONG)
CHECK_TYPE_SIZE("size_t" SIZEOF_SIZE_T)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)

if(NOT (DEFINED NO_SOCKETS))
  check_function_exists(socket SOCKETS_FOUND)
  if(NOT SOCKETS_FOUND)
    set(NO_SOCKETS TRUE)
  endif()
endif()
if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(WINDOWS TRUE)
endif()

if(NO_SOCKETS AND !WINDOWS)
  set(FM_SOCKET_IO null_io.c)
elseif(WINDOWS)
  set(FM_SOCKET_IO "server_acts.c;nt_io.c")
  set(FFS_FILE_IO nt_io.c)
else()
  set(FM_SOCKET_IO "server_acts.c;unix_io.c")
  set(FFS_FILE_IO unix_io.c)
endif()

set(FM_SRC_LIST
   fm_formats.c fm_dump.c lookup3.c string_conversion.c fm_get.c xml.c
   ${FM_SOCKET_IO} ${FFS_FILE_IO})
foreach(FM_SRC ${FM_SRC_LIST})
  list(APPEND FM_MASTER_SRC_LIST fm/${FM_SRC} )
endforeach()

set(COD_SRC_LIST cg.c standard.c)
foreach(COD_SRC ${COD_SRC_LIST})
  list(APPEND COD_MASTER_SRC_LIST cod/${COD_SRC})
endforeach()

set(FFS_SRC_LIST
  ffs.c ffs_formats.c ffs_conv.c ffs_gen.c ffs_file.c evol.c
  ffs_marshal.c)
foreach(FFS_SRC ${FFS_SRC_LIST})
  list(APPEND FFS_MASTER_SRC_LIST ffs/${FFS_SRC})
endforeach()
list(APPEND FFS_MASTER_SRC_LIST version.c)

find_package(BISON)
find_package(FLEX)

if(NOT BISON_FOUND OR NOT FLEX_FOUND)
  include(BisonFlexSub)
  SETUP_BISON_FLEX_SUB()
else()
  BISON_TARGET(CODParser
    cod/cod.y
    ${CMAKE_CURRENT_BINARY_DIR}/cod.tab.c
    COMPILE_FLAGS -d)
    set(HAVE_COD_PARSER_H)
  FLEX_TARGET(CODScanner
    cod/cod.l
    ${CMAKE_CURRENT_BINARY_DIR}/lex.yy.c)
  ADD_FLEX_BISON_DEPENDENCY(CODScanner CODParser)
endif()

find_package(Perl REQUIRED)

add_custom_command(
  OUTPUT "cod_node.c" 
  COMMAND
    ${PERL_EXECUTABLE} 
      ${CMAKE_CURRENT_SOURCE_DIR}/cod/struct.pl
      ${CMAKE_CURRENT_SOURCE_DIR}/cod/cod.structs
  DEPENDS
    lex.yy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/cod/cod.structs
    ${CMAKE_CURRENT_SOURCE_DIR}/cod/struct.pl)

add_custom_target(docs ALL DEPENDS "cod_node.c")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
   SET_SOURCE_FILES_PROPERTIES(cod.tab.c PROPERTIES COMPILE_FLAGS -Wno-unused-function)
endif()

add_custom_target(generated DEPENDS cod_node.c ${BISON_CODParser_OUTPUT_SOURCE})

list(APPEND COD_MASTER_SRC_LIST
  ${BISON_CODParser_OUTPUT_SOURCE}
  ${CMAKE_CURRENT_BINARY_DIR}/cod_node.c)
  
list(APPEND FFS_MASTER_SRC_LIST ${FM_MASTER_SRC_LIST} ${COD_MASTER_SRC_LIST})

set(FFS_LIBRARY_PREFIX "" CACHE STRING
  "Prefix to prepend to the output library name")
mark_as_advanced(FFS_LIBRARY_PREFIX)

add_library(ffs ${FFS_MASTER_SRC_LIST})
set_target_properties(ffs PROPERTIES
  OUTPUT_NAME ${FFS_LIBRARY_PREFIX}ffs
  VERSION ${FFS_VERSION}
  SOVERSION ${FFS_VERSION_MAJOR})
add_library(ffs::ffs ALIAS ffs)

if(UNIX)
  target_link_libraries(ffs PRIVATE m)
  list(APPEND _pkg_config_private_libs m)
endif()
if(SHARED_LIBS_SUPPORTED)
  unset(NO_DYNAMIC_LINKING)
  target_link_libraries(ffs PRIVATE ${CMAKE_DL_LIBS})
  list(APPEND _pkg_config_private_libs ${CMAKE_DL_LIBS})
else()
  set(NO_DYNAMIC_LINKING TRUE)
endif()
#if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(NO_DYNAMIC_LINKING TRUE)
#endif()

target_include_directories(ffs PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fm>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cod>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ffs>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/ffs>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

find_package(dill 3.0.0)

option(FFS_USE_DILL "Enable Dill code generation" ${DILL_FOUND})

if (DILL_FOUND)
  set(HAVE_DILL_H ON)
  option(FFS_CONVERSION_GENERATION_DEFAULT "Enable dynamic code generation" ${DILL_HAS_NATIVE_DCG})
  mark_as_advanced(FFS_CONVERSION_GENERATION_DEFAULT)
  set(DO_DCG TRUE)
  target_link_libraries(ffs PRIVATE dill::dill)
  list(APPEND _pkg_config_private_reqs "dill >= ${dill_VERSION}")
  set(COD_INSTALL_FILE "cod/cod.h")
endif()

option(FFS_USE_ATL "Enable the use of ATL" ON)
if(FFS_USE_ATL)
  find_package(atl 2.2.1 REQUIRED)
  set(HAVE_ATL_H ON)
  target_link_libraries(ffs PUBLIC atl::atl)
  list(APPEND _pkg_config_private_reqs "atl >= ${atl_VERSION}")
endif()

CHECK_INCLUDE_FILE(malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILE(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILE(netdb.h HAVE_NETDB_H)
CHECK_INCLUDE_FILE(netinet/in.h HAVE_NETINET_IN_H)
CHECK_INCLUDE_FILE(sockLib.h HAVE_SOCKLIB_H)
CHECK_INCLUDE_FILE(stdarg.h STDC_HEADERS)
CHECK_INCLUDE_FILE(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILE(string.h HAVE_STRING_H)
CHECK_INCLUDE_FILE(sys/select.h HAVE_SYS_SELECT_H)
CHECK_INCLUDE_FILE(sys/socket.h HAVE_SYS_SOCKET_H)
CHECK_INCLUDE_FILE(sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE(sys/times.h HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILE(sys/uio.h HAVE_SYS_UIO_H)
CHECK_INCLUDE_FILE(sys/un.h HAVE_SYS_UN_H)
CHECK_INCLUDE_FILE(unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILE(windows.h HAVE_WINDOWS_H)
CHECK_INCLUDE_FILE(winsock.h HAVE_WINSOCK_H)

if(SIZEOF_SIZE_T EQUAL SIZEOF_INT)
  set(UIO_SIZE_T_TYPE "unsigned int")
else()
  set(UIO_SIZE_T_TYPE "size_t")
endif()

configure_file(ffs/ffs.h.in ffs/ffs.h @ONLY)
configure_file(CTestCustom.ctest.in CTestCustom.ctest @ONLY)

if(SIZEOF_LONG EQUAL 8)
  set(DATA_LEN_TYPE "long")
elseif(SIZEOF_LONG_LONG EQUAL 8)
  set(DATA_LEN_TYPE "long long")
else()
  set(DATA_LEN_TYPE "undefined")
endif()

CHECK_FUNCTION_EXISTS(strtof HAVE_STRTOF)
CHECK_FUNCTION_EXISTS(strtod HAVE_STRTOD)
CHECK_FUNCTION_EXISTS(strtold HAVE_STRTOLD)
CHECK_FUNCTION_EXISTS(getdomainname HAVE_GETDOMAINNAME)

set(CMAKE_EXTRA_INCLUDE_FILES sys/socket.h)
CHECK_TYPE_SIZE("struct iovec" IOVEC_DEFINE)
unset(CMAKE_EXTRA_INCLUDE_FILES)
CHECK_STRUCT_HAS_MEMBER("struct iovec" iov_base sys/uio.h HAS_IOV_BASE_IOVEC)

if(NOT HAS_IOV_BASE_IOVEC)
  set(NEED_IOVEC_DEFINE TRUE)
endif()

if(NOT DEFINED FORMAT_SERVER_HOST)
  set(FORMAT_SERVER_HOST "formathost.cercs.gatech.edu")
endif()

if(NOT DEFINED FORMAT_SERVICE_DOMAIN)
  set(FORMAT_SERVICE_DOMAIN "")
endif()

# Setup pkgconfig
option(FFS_INSTALL_PKGCONFIG "Install FFS pkgconfig files" ON)
mark_as_advanced(FFS_INSTALL_PKGCONFIG)
if(FFS_INSTALL_PKGCONFIG)
  set(__pkg_config_private_libs ${_pkg_config_private_libs})
  set(_pkg_config_private_libs)
  foreach(L ${__pkg_config_private_libs})
    if(L MATCHES "(.*)/?lib(.*)\\.")
      if(CMAKE_MATCH_1)
        list(APPEND _pkg_config_private_libs "-L${CMAKE_MATCH_1}")
      endif()
      list(APPEND _pkg_config_private_libs "-l${CMAKE_MATCH_2}")
    elseif(L MATCHES "^-")
      list(APPEND _pkg_config_private_libs "${L}")
    else()
      list(APPEND _pkg_config_private_libs "-l${L}")
    endif()
  endforeach()
  string(REPLACE ";" " " _pkg_config_private_libs "${_pkg_config_private_libs}")
  string(REPLACE ";" " " _pkg_config_private_reqs "${_pkg_config_private_reqs}")
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/ffs.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/ffs.pc
    @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ffs.pc
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
    COMPONENT ${FFS_HEADER_COMPONENT})
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/ffs-config.in
    ${CMAKE_CURRENT_BINARY_DIR}/ffs-config
    @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ffs-config
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    COMPONENT ${FFS_HEADER_COMPONENT})
endif()

option(FFS_INSTALL_HEADERS "Install FFS header files" ON)
mark_as_advanced(FFS_INSTALL_HEADERS)
if(FFS_INSTALL_HEADERS)
  install(FILES fm/fm.h ${COD_INSTALL_FILE} ${CMAKE_CURRENT_BINARY_DIR}/ffs/ffs.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT ${FFS_HEADER_COMPONENT})
endif()

set(namelink_component_args)
if(NOT CMAKE_VERSION VERSION_LESS 3.12)
  set(namelink_component_args NAMELINK_COMPONENT ${FFS_HEADER_COMPONENT})
endif()
install(TARGETS ffs EXPORT ffs-targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT ${FFS_RUNTIME_COMPONENT}
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT ${FFS_LIBRARY_COMPONENT} ${namelink_component_args}
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT ${FFS_ARCHIVE_COMPONENT}
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ffs" COMPONENT ${FFS_HEADER_COMPONENT})

if(CMAKE_C_COMPILER_ID MATCHES Intel)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -shared-intel")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# Add all targets to the build-tree export set
export(TARGETS ffs NAMESPACE ffs::
  FILE "${PROJECT_BINARY_DIR}/ffs-targets.cmake")

# Create the ffs-config.cmake and ffs-config-version files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/ffs-config-version.cmake"
  COMPATIBILITY SameMajorVersion)
configure_file(ffs-config.cmake.in
  "${PROJECT_BINARY_DIR}/ffs-config.cmake" @ONLY)

# Install the ffs-config.cmake and ffs-config-version.cmake
install(
  FILES
    "${PROJECT_BINARY_DIR}/ffs-config.cmake"
    "${PROJECT_BINARY_DIR}/ffs-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_CMAKEDIR}" COMPONENT ${FFS_HEADER_COMPONENT})

# Install the export set for use with the install-tree
install(EXPORT ffs-targets NAMESPACE ffs::
  DESTINATION "${CMAKE_INSTALL_CMAKEDIR}" COMPONENT ${FFS_HEADER_COMPONENT})

include(CTest)
if(BUILD_TESTING)
  enable_testing()
endif()

# display status message for important variables
option(FFS_LIBRARIES_ONLY "Whether or not to build the associated executables"
  OFF)
option(FFS_QUIET "Suppress summary output" OFF)
if(NOT FFS_QUIET)
  message(STATUS)
  message(STATUS "-----------------------------------------------------------------------------")
  message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
  message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
  message(STATUS "FFS_USE_DILL = ${FFS_USE_DILL}")
  message(STATUS "FFS_USE_ATL = ${FFS_USE_ATL}")
  message(STATUS "FFS_LIBRARIES_ONLY = ${FFS_LIBRARIES_ONLY}")
  message(STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}")
  message(STATUS "BUILD_TESTING = ${BUILD_TESTING}")
  message(STATUS "Change a value with: cmake -D<Variable>=<Value>")
  message(STATUS "-----------------------------------------------------------------------------")
endif()

add_subdirectory(fm)
add_subdirectory(cod)
add_subdirectory(ffs)
