cmake_minimum_required(VERSION 3.5)

# The directory label is used for CDash to treat EVPath as a subproject of
# GTKorvo
set(CMAKE_DIRECTORY_LABELS EVPath)

project(EVPath VERSION 4.5.0 LANGUAGES C CXX)

# Enable <PackageName>_ROOT variables for dependency searching
# CMake v3.12
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# Enable imported targets as library dependencies for CHECK_INCLUDE_FILE
if(POLICY CMP0075)
  cmake_policy(SET CMP0075 NEW)
endif()

# Let option(...) setings work with non-cache variables.  Useful for including
# EVPath as a third party dependency in another project.
# CMake v3.13
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif()

# Some boilerplate to setup nice output directories
include(GNUInstallDirs)
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/EVPath
  CACHE STRING "Installation CMake subdirectory")
mark_as_advanced(CMAKE_INSTALL_CMAKEDIR)

set(EVPATH_INSTALL_MODULE_DIR "${CMAKE_INSTALL_LIBDIR}"
  CACHE STRING "Directory in which to install EVPath transport modules")
mark_as_advanced(EVPATH_INSTALL_MODULE_DIR)

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
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT DEFINED EVPATH_RUNTIME_COMPONENT)
  set(EVPATH_RUNTIME_COMPONENT bin)
endif()
if(NOT DEFINED EVPATH_LIBRARY_COMPONENT)
  set(EVPATH_LIBRARY_COMPONENT shlib)
endif()
if(NOT DEFINED EVPATH_ARCHIVE_COMPONENT)
  set(EVPATH_ARCHIVE_COMPONENT lib)
endif()
if(NOT DEFINED EVPATH_HEADER_COMPONENT)
  set(EVPATH_HEADER_COMPONENT dev)
endif()
if(NOT DEFINED EVPATH_MODULE_COMPONENT_PREFIX)
  set(EVPATH_MODULE_COMPONENT_PREFIX "")
endif()

# Default to a RelWithDebInfo build if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE RelWithDebInfo)
endif()

if(NOT MSVC)
  set(CMAKE_C_STANDARD 99)
  set(CMAKE_C_STANDARD_REQUIRED True)
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

set(CPACK_DEBIAN_PACKAGE_DEPENDS "dill, atl, ffs")
set(CPACK_RPM_PACKAGE_REQUIRES "dill, atl, ffs")
set(ENABLE_SOMETHING AUTO CACHE STRING "Enable SOMETHING support")  #
set(EVPATH_DEFAULT_PORT_RANGE  "26000:26100" CACHE STRING "Default port range for IP communication")
mark_as_advanced(EVPATH_DEFAULT_PORT_RANGE)

include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CheckTypeSize)
include(CheckStructHasMember)
include(CheckCSourceRuns)
include(TestBigEndian)
include(CheckBrokenTitanCompiler)

include(CTest)
configure_file(CTestCustom.ctest.in CTestCustom.ctest @ONLY)

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

if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall ")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall ")
endif()

add_library(EVPath
  cm.c cm_control.c cm_formats.c cm_util.c cm_transport.c 
  cm_lock.c cm_perf.c cm_pbio.c cm_interface.c version.c
  cm_threadio.c cm_evol.c evp.c response.c metrics.c 
  dlloader.c ip_config.c chr_time.c
  revp.c evp_compat.c thin_server.c evp_threads.c ev_dfg.c)
add_library(EVPath::EVPath ALIAS EVPath)
add_library(evpath ALIAS EVPath)
target_include_directories(EVPath PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

function(_pkg_get_target_prefix tgt out_prefix)
  set(prefix)
  get_property(inc TARGET ${tgt} PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  get_filename_component(prefix "${inc}" PATH)
  set(${out_prefix} "${prefix}" PARENT_SCOPE)
endfunction()

set(_pkg_config_pfxs
  "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/pkgconfig")
set(_pkg_config_reqs)
set(_pkg_config_private_reqs)
set(_pkg_config_libs)
set(_pkg_config_private_libs)

include(CheckCSourceCompiles)

# Avoid recent compilers to optimizing out libm function calls
set(CMAKE_REQUIRED_FLAGS "-O0")
set(LIBM_TEST_SOURCE "#include<math.h>\nfloat f; int main(){sqrt(f);return 0;}")
check_c_source_compiles("${LIBM_TEST_SOURCE}" HAVE_MATH)
if(NOT HAVE_MATH)
  set(CMAKE_REQUIRED_LIBRARIES m)
  check_c_source_compiles("${LIBM_TEST_SOURCE}" HAVE_LIBM_MATH)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(NOT HAVE_LIBM_MATH)
    message(FATAL_ERROR "Unable to use C math library functions")
  endif()
    
  target_link_libraries(EVPath PRIVATE m)
  list(APPEND _pkg_config_private_libs m)
endif()

set(EVPATH_LIBRARY_PREFIX "" CACHE STRING
  "Prefix to prepend to the output library name")
mark_as_advanced(EVPATH_LIBRARY_PREFIX)
set_target_properties(EVPath PROPERTIES
  OUTPUT_NAME ${EVPATH_LIBRARY_PREFIX}evpath)

set(IPCONFIG_ENVVAR_PREFIX "CM_" CACHE STRING
  "Prefix to prepend to the networking environment variable names")
mark_as_advanced(IPCONFIG_ENVVAR_PREFIX)
if (NOT EVPATH_LIBRARY_PREFIX STREQUAL "")
   STRING (TOUPPER ${EVPATH_LIBRARY_PREFIX} IPCONFIG_ENVVAR_PREFIX)
endif()
   
mark_as_advanced(EVPATH_LIBRARY_PREFIX)

set(EVPATH_TRANSPORT_DEP_LIBS)

if(NOT DEFINED CM_SELF_FORMATS)
  set(CM_SELF_FORMATS 1)
endif()

if (DEFINED EVPATH_NO_RDMA) 
  option(EVPATH_NO_RDMA "Don't use any RDMA library"
    "${EVPATH_NO_RDMA}")
  if (EVPATH_NO_RDMA)
    SET (EVPATH_USE_LIBFABRIC OFF)
    SET (EVPATH_USE_NNTI OFF)
    SET (EVPATH_USE_IBVERBS OFF)
  endif()
endif()

cmake_dependent_option(EVPATH_TRANSPORT_MODULES
  "Enable stand-alone transport modules" ${BUILD_SHARED_LIBS}
  "SHARED_LIBS_SUPPORTED" OFF
)

list(INSERT CMAKE_PREFIX_PATH 0 ${CMAKE_INSTALL_PREFIX})
find_package(atl 2.2.1 REQUIRED)
find_package(ffs 3.0.0 REQUIRED)
_pkg_get_target_prefix(atl::atl atl_PREFIX)
_pkg_get_target_prefix(ffs::ffs ffs_PREFIX)
list(APPEND _pkg_config_pfxs "${atl_PREFIX}" "${ffs_PREFIX}")

if(TARGET ffs) # The unaliased target exists so it's built by us
  set(HAVE_COD_H TRUE CACHE INTERNAL "" FORCE)
else()
  set(CMAKE_REQUIRED_LIBRARIES ffs::ffs)
  CHECK_INCLUDE_FILE(cod.h HAVE_COD_H)
  unset(CMAKE_REQUIRED_LIBRARIES)
endif()
if (HAVE_COD_H)
   find_package(dill 3.0.0 REQUIRED)
   target_link_libraries(EVPath PUBLIC atl::atl ffs::ffs PRIVATE dill::dill)
   _pkg_get_target_prefix(dill::dill dill_PREFIX)
   list(APPEND _pkg_config_private_reqs "dill >= ${dill_VERSION}")
   list(APPEND _pkg_config_pfxs "${dill_PREFIX}")
else()
   target_link_libraries(EVPath PUBLIC atl::atl ffs::ffs)
endif()

find_package(Threads)
if(Threads_FOUND)
  target_link_libraries(EVPath PRIVATE ${CMAKE_THREAD_LIBS_INIT})
  list(APPEND _pkg_config_private_libs ${CMAKE_THREAD_LIBS_INIT})
endif()

if(DEFINED EVPATH_USE_NVML)
  option(EVPATH_USE_NVML "Use the NVidia management library"
    "${EVPATH_USE_NVML}")
  if(EVPATH_USE_NVML)
    find_package(NVML REQUIRED)
  endif()
else()
  find_package(NVML)
  option(EVPATH_USE_NVML "Use the NVidia management library" ${NVML_FOUND})
endif()
if(NVML_FOUND)
  if(NOT BUILD_SHARED_LIBS)
    list(APPEND EXTRA_INSTALL_PACKAGES NVML)
  endif()
  target_link_libraries(EVPath PRIVATE nvml::nvml)
  list(APPEND _pkg_config_private_libs ${NVML_LIBRARIES})
endif()

include(CheckSymbolExists)
CHECK_SYMBOL_EXISTS(clock_gettime time.h HAVE_CLOCK_GETTIME)

if(NOT HAVE_CLOCK_GETTIME)
  include(CheckLibraryExists)
  CHECK_LIBRARY_EXISTS(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME)
  if(HAVE_CLOCK_GETTIME)
    target_link_libraries(EVPath PRIVATE rt)
    list(APPEND _pkg_config_private_libs rt)
  endif()
endif()

CHECK_INCLUDE_FILE(sys/epoll.h HAVE_SYS_EPOLL_H)
set(CM_DEFAULT_SELECT "select")
if (HAVE_SYS_EPOLL_H)
   set(CM_DEFAULT_SELECT "epoll")
endif()
set(CM_DEFAULT_TRANSPORT "sockets")

# This is a dummy interace library to allow the transport modules to depend on
# the EVPath headers without causing an unnecessary link dependency
add_library(evpath_headers INTERFACE)
set_target_properties(evpath_headers PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES
    "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR};$<TARGET_PROPERTY:ffs::ffs,INTERFACE_INCLUDE_DIRECTORIES>;$<TARGET_PROPERTY:atl::atl,INTERACE_INCLUDE_DIRECTORIES>"
)

set(EVPATH_TRANSPORT_TARGETS)
if(EVPATH_TRANSPORT_MODULES)
  target_link_libraries(EVPath PRIVATE ${CMAKE_DL_LIBS})

  add_library(cmselect MODULE cmselect.c)
  add_library(cmsockets MODULE cmsockets.c ip_config.c)
  add_library(cmudp MODULE cmudp.c)
  add_library(cmmulticast MODULE cmmulticast.c)

  foreach(M cmselect cmsockets cmudp cmmulticast)
    set_target_properties(${M} PROPERTIES
      OUTPUT_NAME ${EVPATH_LIBRARY_PREFIX}${M}
      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${EVPATH_INSTALL_MODULE_DIR}
      WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    install(TARGETS ${M} DESTINATION ${EVPATH_INSTALL_MODULE_DIR} COMPONENT ${EVPATH_MODULE_COMPONENT_PREFIX}${M})
  endforeach()

  target_link_libraries(cmselect PRIVATE evpath_headers atl::atl)
  target_link_libraries(cmsockets PRIVATE evpath_headers atl::atl)
  target_link_libraries(cmudp PRIVATE evpath_headers atl::atl)
  target_link_libraries(cmmulticast PRIVATE evpath_headers atl::atl)

  list(APPEND EVPATH_TRANSPORT_TARGETS cmselect cmsockets cmudp cmmulticast)

  if (HAVE_SYS_EPOLL_H)
     add_library(cmepoll MODULE cmepoll.c)
     set_target_properties(cmepoll PROPERTIES
          OUTPUT_NAME ${EVPATH_LIBRARY_PREFIX}cmepoll
          LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${EVPATH_INSTALL_MODULE_DIR})
     target_link_libraries(cmepoll PRIVATE evpath_headers atl::atl)
     list(APPEND EVPATH_TRANSPORT_TARGETS cmepoll)
     install(TARGETS cmepoll DESTINATION ${EVPATH_INSTALL_MODULE_DIR} COMPONENT ${EVPATH_MODULE_COMPONENT_PREFIX}cmepoll)
  endif()
else()
  # dummy dyn linking in cm.c and cm_transport.c
  set(NO_DYNAMIC_LINKING 1)

  target_sources(EVPath PRIVATE cmsockets.c)
  target_sources(EVPath PRIVATE cmselect.c)
  if (HAVE_SYS_EPOLL_H)
    target_sources(EVPath PRIVATE cmepoll.c)
  endif()
  target_sources(EVPath PRIVATE cmudp.c)
  target_sources(EVPath PRIVATE cmmulticast.c)
endif()

set(RUN_NNTI_TESTS FALSE)

# enet transport
set(RUN_ENET_TESTS FALSE)
if(DEFINED EVPATH_USE_ENET)
  option(EVPATH_USE_ENET "Build the enet transport" "${EVPATH_USE_ENET}")
  if(EVPATH_USE_ENET)
    find_package(enet 1.3.13 REQUIRED)
  endif()
else()
  find_package(enet 1.3.13 QUIET)
  option(EVPATH_USE_ENET "Build the enet transport" ${ENET_FOUND})
endif()
if(ENET_FOUND)
  set(RUN_ENET_TESTS TRUE)
  if(EVPATH_TRANSPORT_MODULES)
    add_library(cmenet MODULE cmenet.c ip_config.c)
    set_target_properties(cmenet PROPERTIES
      OUTPUT_NAME ${EVPATH_LIBRARY_PREFIX}cmenet
      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${EVPATH_INSTALL_MODULE_DIR})

    target_link_libraries(cmenet PRIVATE evpath_headers atl::atl enet::enet)
    list(APPEND EVPATH_TRANSPORT_TARGETS cmenet)
    install(TARGETS cmenet DESTINATION ${EVPATH_INSTALL_MODULE_DIR} COMPONENT ${EVPATH_MODULE_COMPONENT_PREFIX}cmenet)
  else()
    target_sources(EVPath PRIVATE cmenet.c)
    target_link_libraries(EVPath PRIVATE enet::enet)
    _pkg_get_target_prefix(enet::enet enet_PREFIX)
    list(APPEND _pkg_config_pfxs "${enet_PREFIX}")
    list(APPEND _pkg_config_private_reqs "enet >= ${enet_VERSION}")
  endif()
else()
  message(STATUS " - Enet library was not found.  This is not a fatal error, just that the Enet transport will not be built.")
endif()

if (MSVC)
    set(EVPATH_USE_ZPL_ENET FALSE)
endif()
if(NOT (DEFINED EVPATH_USE_ZPL_ENET))
  option(EVPATH_USE_ZPL_ENET "Build the zplenet transport" "OFF")
endif()
if(EVPATH_USE_ZPL_ENET)
  set(RUN_ZPL_ENET_TESTS TRUE)
  set(ZPL_ENET_AVAILABLE TRUE)
  CONFIGURE_FILE(cmenet.c cmzplenet.c COPYONLY)
  SET_SOURCE_FILES_PROPERTIES(cmzplenet.c PROPERTIES COMPILE_FLAGS "-DUSE_ZPL_ENET=1 -DENET_API=static")
  if(EVPATH_TRANSPORT_MODULES)
    add_library(cmzplenet MODULE cmzplenet.c ip_config.c)
    set_target_properties(cmzplenet PROPERTIES
      OUTPUT_NAME ${EVPATH_LIBRARY_PREFIX}cmzplenet
      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${EVPATH_INSTALL_MODULE_DIR})

    target_link_libraries(cmzplenet PRIVATE evpath_headers atl::atl)
    list(APPEND EVPATH_TRANSPORT_TARGETS cmzplenet)
    install(TARGETS cmzplenet DESTINATION ${EVPATH_INSTALL_MODULE_DIR} COMPONENT ${EVPATH_MODULE_COMPONENT_PREFIX}cmzplenet)
  else()
    target_sources(EVPath PRIVATE cmzplenet.c)
  endif()
endif()

# udt4 transport
set(RUN_UDT4_TESTS FALSE)
if(NOT(DEFINED EVPATH_USE_UDT4))
  #UDT4 not particularly functional, default to not building
  set(EVPATH_USE_UDT4 FALSE)
endif()
if(DEFINED EVPATH_USE_UDT4)
  option(EVPATH_USE_UDT4 "Build the udt4 transport" "${EVPATH_USE_UDT4}")
  if(EVPATH_USE_UDT4)
    find_package(udt4 4.11 REQUIRED)
  endif()
else()
  find_package(udt4 4.11)
  option(EVPATH_USE_UDT4 "Build the udt4 transport" ${UDT4_FOUND})
  if (NOT UDT4_FOUND)
    message(STATUS " - Udt4 library was not found.  This is not a fatal error, just that the Udt4 transport will not be built.")
  endif()
endif()
if(UDT4_FOUND)
  set(RUN_UDT4_TESTS FALSE)   # DON'T TEST.  UDT4 NOT USEFUL
  if(EVPATH_TRANSPORT_MODULES)
    add_library(cmudt4 MODULE cmudt4.cpp ip_config.c)
    set_target_properties(cmudt4 PROPERTIES
      OUTPUT_NAME ${EVPATH_LIBRARY_PREFIX}cmudt4
      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${EVPATH_INSTALL_MODULE_DIR})

    target_link_libraries(cmudt4 PRIVATE evpath_headers atl::atl udt4::udt4)
    list(APPEND EVPATH_TRANSPORT_TARGETS cmudt4)
    install(TARGETS cmudt4 DESTINATION ${EVPATH_INSTALL_MODULE_DIR} COMPONENT ${EVPATH_MODULE_COMPONENT_PREFIX}cmudt4)
  else()
    target_sources(EVPath PRIVATE cmudt4.cpp)
    target_link_libraries(EVPath PRIVATE udt4::udt4)
    _pkg_get_target_prefix(udt4::udt4 udt4_PREFIX)
    list(APPEND _pkg_config_pfxs "${udt4_PREFIX}")
    list(APPEND _pkg_config_private_reqs "udt4 >= ${udt4_VERSION}")
  endif()
endif()

# libfabric transport
if(DEFINED EVPATH_USE_LIBFABRIC)
  option(EVPATH_USE_LIBFABRIC "Build the libfabric transport"
    "${EVPATH_USE_LIBFABRIC}")
  if(EVPATH_USE_LIBFABRIC)
    find_package(LIBFABRIC REQUIRED)
  endif()
else()
  find_package(LIBFABRIC)
  option(EVPATH_USE_LIBFABRIC "Build the libfabric transport"
    ${LIBFABRIC_FOUND})
endif()
if(LIBFABRIC_FOUND AND NOT EVPATH_NO_RDMA) 
  if(EVPATH_TRANSPORT_MODULES)
    add_library(cmfabric MODULE cmfabric.c ip_config.c)
    set_target_properties(cmfabric PROPERTIES
      OUTPUT_NAME ${EVPATH_LIBRARY_PREFIX}cmfabric
      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${EVPATH_INSTALL_MODULE_DIR})
    target_link_libraries(cmfabric PRIVATE
      evpath_headers atl::atl libfabric::libfabric)
    list(APPEND EVPATH_TRANSPORT_TARGETS cmfabric)
    install(TARGETS cmfabric DESTINATION ${EVPATH_INSTALL_MODULE_DIR} COMPONENT ${EVPATH_MODULE_COMPONENT_PREFIX}cmfabric)
  else()
    if(NOT BUILD_SHARED_LIBS)
      list(APPEND EXTRA_INSTALL_PACKAGES LIBFABRIC)
    endif()
    target_sources(EVPath PRIVATE cmfabric.c)
    target_link_libraries(EVPath PRIVATE libfabric::libfabric)
    list(APPEND _pkg_config_pfxs ${LIBFABRIC_LIBDIR}/pkgconfig)
    list(APPEND _pkg_config_private_reqs "libfabric >= ${LIBFABRIC_VERSION}")
  endif()
else()
  if (EVPATH_USE_LIBFABRIC)
    message(STATUS " - LibFabric package was not found.  This is not a fatal error, just that the fabric transport will not be built.")
  endif()
endif()

# ibverbs transport
set(RUN_IB_TESTS FALSE)
if(DEFINED EVPATH_USE_IBVERBS)
  option(EVPATH_USE_IBVERBS "Build the libfabric transport"
    "${EVPATH_USE_IBVERBS}")
  if(EVPATH_USE_IBVERBS)
    find_package(IBVERBS REQUIRED)
  endif()
else()
  find_package(IBVERBS)
  option(EVPATH_USE_IBVERBS "Build the libfabric transport" ${IBVERBS_FOUND})
endif()
set(HAVE_IBVERBS ${IBVERBS_FOUND})
if(IBVERBS_FOUND AND NOT EVPATH_NO_RDMA)
  if(BUILD_TESTING)
    if(NOT CMAKE_CROSSCOMPILING)
      message(STATUS "Check MEMLOCK rlimit for IB tests")
      file(READ "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GOOD_MEMLOCK_LIMIT.c" _SOURCE)
      CHECK_C_SOURCE_RUNS("${_SOURCE}" GOOD_MEMLOCK_LIMIT)
      if(GOOD_MEMLOCK_LIMIT)
        message(STATUS "Check MEMLOCK rlimit for IB tests - yes" )
        set(RUN_IB_TESTS TRUE)
        set(RUN_NNTI_TESTS FALSE)
      else()
        message(STATUS "Check MEMLOCK rlimit for IB tests - no")
	message(STATUS "  Building InfiniBand transport, but current")
        message(STATUS "  RLIMIT_MEMLOCK value will prevent successful runs.")
        message(STATUS "  Not running IB transport tests.")
      endif()
    else()
      set(RUN_IB_TESTS TRUE)
      set(RUN_NNTI_TESTS FALSE)
    endif()
  endif()

  if(WIN32)
  target_link_libraries(EVPath wsock32 ws2_32)
  endif()

  if(EVPATH_TRANSPORT_MODULES)
    add_library(cmib MODULE cmib.c)
    set_target_properties(cmib PROPERTIES
      OUTPUT_NAME ${EVPATH_LIBRARY_PREFIX}cmib
      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${EVPATH_INSTALL_MODULE_DIR})
    target_link_libraries(cmib PRIVATE
      evpath_headers atl::atl ${IBVERBS_LIBRARIES})
    list(APPEND EVPATH_TRANSPORT_TARGETS cmib)
    install(TARGETS cmib DESTINATION ${EVPATH_INSTALL_MODULE_DIR} COMPONENT ${EVPATH_MODULE_COMPONENT_PREFIX}cmib)
  else()
    target_sources(EVPath PRIVATE cmib.c)
    target_link_libraries(EVPath PRIVATE ${IBVERBS_LIBRARIES})
    list(APPEND _pkg_config_private_libs ${IBVERBS_LIBRARIES})
  endif()
  set(IB_FOUND 1)
endif()

if(NOT (DEFINED CercsArch))
  execute_process(
    COMMAND cercs_arch
    OUTPUT_VARIABLE CercsArch
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  mark_as_advanced(CercsArch)
endif()

if(IBVERBS_FOUND)
  if(DEFINED EVPATH_USE_NNTI)
    option(EVPATH_USE_NNTI "Build with nnti support" "${EVPATH_USE_NNTI}")
    if(EVPATH_USE_NNTI)
      find_package(NNTI REQUIRED)
    endif()
  else()
    find_package(NNTI)
    option(EVPATH_USE_NNTI "Build the nnti transport" ${NNTI_FOUND})
  endif()
  if(NNTI_FOUND AND NOT EVPATH_NO_RDMA)
    if(CercsArch STREQUAL "ppc64")
      set(RUN_NNTI_TESTS FALSE)
    endif()
    if(EVPATH_TRANSPORT_MODULES)
      add_library(cmnnti MODULE cmnnti.c)
      set_target_properties(cmnnti PROPERTIES
        OUTPUT_NAME ${EVPATH_LIBRARY_PREFIX}cmnnti
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${EVPATH_INSTALL_MODULE_DIR})
      target_link_libraries(cmnnti PRIVATE
        evpath_headers atl::atl nnti::nnti enet::enet)
      list(APPEND EVPATH_TRANSPORT_TARGETS cmnnti)
    install(TARGETS cmnnti DESTINATION ${EVPATH_INSTALL_MODULE_DIR} COMPONENT ${EVPATH_MODULE_COMPONENT_PREFIX}cmnnti)
    else()
      if(NOT BUILD_SHARED_LIBS)
        list(APPEND EXTRA_INSTALL_PACKAGES NNTI)
      endif()
      target_sources(EVPath PRIVATE cmnnti.c)
      target_link_libraries(EVPath PRIVATE nnti::nnti enet::enet atl::atl)
      list(APPEND _pkg_config_private_libs ${NNTI_LIBRARIES})
    endif()
  endif()
else()
  set(EVPATH_USE_NNTI FALSE CACHE INTERNAL "Build the nnti transport" FORCE)
endif()

# Install extra find module dependencies
install(DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR} COMPONENT ${EVPATH_HEADER_COMPONENT}
  FILES_MATCHING PATTERN "Find*.cmake" PATTERN "CMake*.cmake"
)

add_custom_command(
  OUTPUT "cm_interface.c" "revp.c" "revpath.h"
  SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/evpath.h ${CMAKE_CURRENT_SOURCE_DIR}/ev_dfg.h
  COMMAND perl
    ${CMAKE_CURRENT_SOURCE_DIR}/gen_interface.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/evpath.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ev_dfg.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cm_schedule.h
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/gen_interface.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/evpath.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ev_dfg.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cm_schedule.h
)

if(Threads_FOUND AND CMAKE_USE_PTHREADS_INIT)
  set(USE_PTHREADS TRUE)
endif()

CHECK_INCLUDE_FILE(hostlib.h HAVE_HOSTLIB_H)
CHECK_INCLUDE_FILE(malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILE(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILE(netdb.h HAVE_NETDB_H)
CHECK_INCLUDE_FILE(sockLib.h HAVE_SOCKLIB_H)
CHECK_INCLUDE_FILE(stdarg.h STDC_HEADERS)
CHECK_INCLUDE_FILE(stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILE(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILE(string.h HAVE_STRING_H)
CHECK_INCLUDE_FILE(sys/select.h HAVE_SYS_SELECT_H)
CHECK_INCLUDE_FILE(sys/socket.h HAVE_SYS_SOCKET_H)
CHECK_INCLUDE_FILE(sys/sockio.h HAVE_SYS_SOCKIO_H)
CHECK_INCLUDE_FILE(sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE(sys/times.h HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILE(sys/uio.h HAVE_SYS_UIO_H)
CHECK_INCLUDE_FILE(sys/un.h HAVE_SYS_UN_H)
CHECK_INCLUDE_FILE(unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILE(windows.h HAVE_WINDOWS_H)
CHECK_INCLUDE_FILE(winsock2.h HAVE_WINSOCK2_H)

CHECK_STRUCT_HAS_MEMBER("struct fd_set" "fds_bits" "sys/select.h" HAVE_FDS_BITS)

CHECK_TYPE_SIZE("int"   SIZEOF_INT)
CHECK_TYPE_SIZE("long"   SIZEOF_LONG)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)

CHECK_FUNCTION_EXISTS(writev HAVE_WRITEV)
CHECK_FUNCTION_EXISTS(uname HAVE_UNAME)
CHECK_FUNCTION_EXISTS(getdomainname HAVE_GETDOMAINNAME)
CHECK_FUNCTION_EXISTS(getloadavg HAVE_GETLOADAVG)
CHECK_FUNCTION_EXISTS(gettimeofday HAVE_GETTIMEOFDAY)
CHECK_FUNCTION_EXISTS(getifaddrs HAVE_GETIFADDRS)

try_compile(HAVE_MAC_SYSCTL
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_mac_sysctl.c
)
try_compile(HAVE_SYSINFO
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_sysinfo.c
)
try_compile(HAVE_SYSCONF
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_sysconf.c
)

if(CMAKE_C_COMPILER_ID MATCHES "Intel") 
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -shared-intel")
endif()

function(_pkg_expand_libs in_libs out_ldflags)
  set(ldflags)
  foreach(L ${in_libs})
    if(L MATCHES "(.*)/?lib(.*).")
      if(CMAKE_MATCH_1)
        list(APPEND ldflags "-L${CMAKE_MATCH_1}")
      endif()
      list(APPEND ldflags "-l${CMAKE_MATCH_2}")
    elseif(L MATCHES "^-")
      list(APPEND ldflags "${L}")
    else()
      list(APPEND ldflags "-l${L}")
    endif()
  endforeach()
  string(REPLACE ";" " " ldflags "${ldflags}")
  set(${out_ldflags} "${ldflags}" PARENT_SCOPE)
endfunction()
string(REPLACE ";" ":" _pkg_config_pfxs "${_pkg_config_pfxs}")
string(REPLACE ";" " " _pkg_config_reqs "${_pkg_config_reqs}")
string(REPLACE ";" " " _pkg_config_private_reqs "${_pkg_config_private_reqs}")
_pkg_expand_libs("${_pkg_config_libs}" _pkg_config_libs)
_pkg_expand_libs("${_pkg_config_private_libs}" _pkg_config_private_libs)

# Setup pkgconfig
option(EVPATH_INSTALL_PKGCONFIG "Install EVPath pkgconfig files" ON)
mark_as_advanced(EVPATH_INSTALL_PKGCONFIG)
if(EVPATH_INSTALL_PKGCONFIG)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/evpath.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/evpath.pc
    @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/evpath.pc
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
    COMPONENT ${EVPATH_HEADER_COMPONENT})
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/evpath_config.in
    ${CMAKE_CURRENT_BINARY_DIR}/evpath_config
    @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/evpath_config
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    COMPONENT ${EVPATH_HEADER_COMPONENT})
endif()

option(EVPATH_INSTALL_HEADERS "Install EVPath header files" ON)
mark_as_advanced(EVPATH_INSTALL_HEADERS)
if(EVPATH_INSTALL_HEADERS)
  install(
    FILES
      evpath.h
      cm_schedule.h
      ev_dfg.h
      cm_transport.h
      ${CMAKE_CURRENT_BINARY_DIR}/revpath.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT ${EVPATH_HEADER_COMPONENT}
  )
endif()

set(namelink_component_args)
if(NOT CMAKE_VERSION VERSION_LESS 3.12)
  set(namelink_component_args NAMELINK_COMPONENT ${EVPATH_HEADER_COMPONENT})
endif()
install(TARGETS EVPath EXPORT EVPathTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT ${EVPATH_RUNTIME_COMPONENT}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${EVPATH_LIBRARY_COMPONENT} ${namelink_component_args}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${EVPATH_ARCHIVE_COMPONENT}
)

FOREACH(tgt ${EVPATH_TRANSPORT_TARGETS})
   TARGET_COMPILE_DEFINITIONS(${tgt} PRIVATE EVPATH_MODULE_BUILD_DIR="$<TARGET_FILE_DIR:${tgt}>")
ENDFOREACH()

set(config_suffix)
get_property(multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(multi_config)
  set(config_suffix "/$<CONFIG>")
endif()
TARGET_COMPILE_DEFINITIONS(EVPath PRIVATE EVPATH_MODULE_BUILD_DIR="${PROJECT_BINARY_DIR}/${EVPATH_INSTALL_MODULE_DIR}${config_suffix}")
set(EVPATH_MODULE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/${EVPATH_INSTALL_MODULE_DIR})

if(CMAKE_C_COMPILER_ID MATCHES "Intel") 
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -shared-intel")
endif()

find_program(SSH ssh)
if(SSH) 
  set(SSH_PATH "${SSH}")
endif() 

if(BUILD_TESTING) 
  add_executable(cmprobe cmprobe.c)
  target_link_libraries(cmprobe evpath atl::atl)

  enable_testing()

  add_subdirectory(tests)
  add_subdirectory(rtests)
  add_subdirectory(mtests)
  add_subdirectory(dfg_tests)

  find_package(MPI)
  if(MPI_FOUND)
    add_subdirectory(mpi_tests)
  endif()

  add_subdirectory(examples)
endif()

export(EXPORT EVPathTargets NAMESPACE EVPath::
  FILE "${PROJECT_BINARY_DIR}/EVPathTargets.cmake")

# This "exports" all targets which have been put into the export set
# "EVPathTargets".  This means that cmake generates a file with the given
# filename, which can later on be loaded by projects using this package.  This
# file contains add_library(EVPath IMPORTED) statements for each target in the
# export set, so when loaded later on cmake will create "imported" library
# targets from these, which can be used in many ways in the same way as a
# normal library target created via a normal add_library().
install(EXPORT EVPathTargets NAMESPACE EVPath::
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT ${EVPATH_HEADER_COMPONENT})

# Create a EVPathConfig.cmake file. <name>Config.cmake files are searched by
# find_package() automatically. We configure that file so that we can put any
# information we want in it, e.g. version numbers, include directories, etc.
configure_file(
  EVPathConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/EVPathConfig.cmake
  @ONLY
)
configure_file(
  EVPathConfigCommon.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/EVPathConfigCommon.cmake
  @ONLY
)

# Additionally, when cmake has found a EVPathConfig.cmake, it can check for a
# EVPathConfigVersion.cmake in the same directory when figuring out the version
# of the package when a version has been specified in the find_package() call,
# e.g. find_package(EVPath 1.0)
configure_file(
  EVPathConfigVersion.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/EVPathConfigVersion.cmake
  @ONLY
)

# Install these two files into the same directory as the generated exports-file.
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/EVPathConfigCommon.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/EVPathConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT ${EVPATH_HEADER_COMPONENT}
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/EVPathConfigInstall.cmake
  RENAME EVPathConfig.cmake
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT ${EVPATH_HEADER_COMPONENT}
)


configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

if(TEST_INSTALL_DIRECTORY) 
  set(EVPATH_TEST_INSTALL_DIR ${TEST_INSTALL_DIRECTORY})
endif()

# display status message for important variables
option(EVPATH_QUIET "Suppress summary output" OFF)
if(NOT EVPATH_QUIET)
  message(STATUS)
  message(STATUS "-----------------------------------------------------------------------------" )
  message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
  message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
  message(STATUS "BUILD_TESTING = ${BUILD_TESTING}  (options are: ON, OFF)")
  message(STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}  (options are: ON, OFF, default OFF)")
  message(STATUS "EVPATH_USE_NVML = ${EVPATH_USE_NVML}")
  message(STATUS "EVPATH_USE_ENET = ${EVPATH_USE_ENET}")
  message(STATUS "EVPATH_USE_ZPL_ENET = ${EVPATH_USE_ZPL_ENET}")
  message(STATUS "EVPATH_USE_UDT4 = ${EVPATH_USE_UDT4}")
  message(STATUS "EVPATH_USE_LIBFABRIC = ${EVPATH_USE_LIBFABRIC}")
  message(STATUS "EVPATH_USE_IBVERBS = ${EVPATH_USE_IBVERBS}")
  message(STATUS "EVPATH_USE_NNTI = ${EVPATH_USE_NNTI}")
  message(STATUS "EVPATH_TRANSPORT_MODULES = ${EVPATH_TRANSPORT_MODULES} (${EVPATH_TRANSPORT_TARGETS})")
  message(STATUS "Change a value with: cmake -D<Variable>=<Value>" )
  message(STATUS "To turn off all RDMA transports (fabric,ib,nnti) in EVPATH set EVPATH_NO_RDMA=1" )
  message(STATUS "-----------------------------------------------------------------------------" )
endif()
