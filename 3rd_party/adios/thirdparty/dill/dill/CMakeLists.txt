cmake_minimum_required(VERSION 3.5)

# The directory label is used for CDash to treat DILL as a subproject of GTKorvo
set(CMAKE_DIRECTORY_LABELS DILL)

project(DILL VERSION 3.0.0 LANGUAGES C CXX)

# Some boilerplate to setup nice output directories
include(GNUInstallDirs)
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/dill
  CACHE STRING "Installation CMake subdirectory")
mark_as_advanced(CMAKE_INSTALL_CMAKEDIR)

list(INSERT CMAKE_PREFIX_PATH 0 ${CMAKE_INSTALL_PREFIX})
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
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

if(NOT DEFINED DILL_RUNTIME_COMPONENT)
  set(DILL_RUNTIME_COMPONENT bin)
endif()
if(NOT DEFINED DILL_LIBRARY_COMPONENT)
  set(DILL_LIBRARY_COMPONENT shlib)
endif()
if(NOT DEFINED DILL_ARCHIVE_COMPONENT)
  set(DILL_ARCHIVE_COMPONENT lib)
endif()
if(NOT DEFINED DILL_HEADER_COMPONENT)
  set(DILL_HEADER_COMPONENT dev)
endif()
add_definitions(-DDILL_SRC)
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

# -O3 is known to crash some flavors of clang 14.
if(CMAKE_C_COMPILER_ID MATCHES "Clang" AND
   CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL 14.0.0)
  foreach(_var_suffix IN ITEMS "" _DEBUG _RELEASE _RELWITHDEBINFO _MINSIZEREL)
    string(REGEX REPLACE "-O(2|3)" "-O" CMAKE_C_FLAGS${_var_suffix}
      "${CMAKE_C_FLAGS${_var_suffix}}")
  endforeach()
endif()

# Setup shared library defaults.  If explicitly specified somehow, then default 
# to that.  Otherwise base the default on whether or not shared libs are even
# supported.
include(CMakeDependentOption)
get_property(SHARED_LIBS_SUPPORTED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
cmake_dependent_option(BUILD_SHARED_LIBS
  "Build shared libraries (so/dylib/dll)." ${SHARED_LIBS_SUPPORTED}
  "SHARED_LIBS_SUPPORTED" OFF
)
mark_as_advanced(BUILD_SHARED_LIBS)

# Default to a RelWithDebInfo build if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE RelWithDebInfo)
endif()

set(SRC_list
  dill_util.c
  dill_pkg.c
  dill_cplus.c
  dill_virtual.c
  virtual.c
  dill.c
  vm.c
  x86_64_disassembler.c
)

set(TARGET_DEP_INC)
set(TARGET_DEP_LIBS)
set(PKG_DEP_PKG)
set(PKG_DEP_LIBS)

include(CheckTypeSize)
check_type_size("void*" CMAKE_SIZEOF_VOID_P)
check_type_size("long" SIZEOF_LONG)
set(NATIVE_CG TRUE)
unset(NATIVE_ARCH)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "i.86|x86_64|AMD64|amd64")
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(NATIVE_ARCH x86_64)
    set(HOST_X86_64 1)
    set(TEST_PERL_FLAGS "-max_arg=4")
    set(ARCH_FILE x86_64)
  else()
    set(NATIVE_ARCH x86)
    set(HOST_X86 1)
    set(TEST_PERL_FLAGS -no_float -max_arg=2)
    set(ARCH_FILE x86)
  endif()
  if(NOT (CMAKE_SIZEOF_VOID_P EQUAL SIZEOF_LONG))
    set(BASE_OPS_ARG "-msvc_long")
  endif()
  if(BUILD_SHARED_LIBS)
    set(BASE_OPS_ARG "${BASE_OPS_ARG} -build_shared")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm.5.*")
  set(NATIVE_ARCH arm5)
  set(HOST_ARM5 1)
  set(TEST_PERL_FLAGS -max_arg=4)
  set(ARCH_FILE arm5)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm.7.")
  set(NATIVE_ARCH arm6)
  set(HOST_ARM7 1)
  set(TEST_PERL_FLAGS -max_arg=4)
  set(ARCH_FILE arm6)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm.*")
  set(NATIVE_ARCH arm6)
  set(HOST_ARM6 1)
  set(TEST_PERL_FLAGS -max_arg=4)
  set(ARCH_FILE arm6)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
  set(NATIVE_ARCH arm8)
  set(HOST_ARM8 1)
  set(ARCH_FILE arm8)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64le")
  set(NATIVE_ARCH ppc64le)
  set(HOST_PPC64LE 1)
  set(ARCH_FILE ppc64le)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "s64fx")
  set(NATIVE_ARCH sparcv9)
  set(HOST_SPARCV9 1)
  set(ARCH_FILE sparc)
else()
  message(STATUS "NO CG for Processor \"${CMAKE_SYSTEM_PROCESSOR}\"")
  set(NATIVE_CG FALSE)
  set(NATIVE_ARCH UNSUPPORTED)
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Darwin|Linux")
  set(USE_MMAP_CODE_SEG 1)
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
  set(USE_VIRTUAL_PROTECT 1)
  set(USE_WINDOWS_CALLS 1)
endif()

include(TestBigEndian)
test_big_endian(WORDS_BIGENDIAN)

set(HARDFP_AVAILABLE 0)
# Determine float ABI of ARM Linux
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^arm" AND
   CMAKE_SYSTEM_NAME STREQUAL "Linux")
  find_program(READELF readelf)
  if(READELF)
    find_library(MATH_LIBRARY m)
    execute_process(
      COMMAND ${READELF} -A ${MATH_LIBRARY}
      OUTPUT_VARIABLE ELF_ARCH
      ERROR_QUIET
    )
    if(ELF_ARCH MATCHES "Tag_ABI_VFP_args: VFP registers")
      set(HARDFP_AVAILABLE 1)
    endif()
  else()
    message(WARNING "readelf not found; float ABI detection disabled")
  endif()
endif()
set(ARM_HARD_FLOAT ${HARDFP_AVAILABLE})

option(DILL_IGNORE_NATIVE  "Build to do emulation, regardless of architecture"
  OFF)
option(DILL_NATIVE_ONLY  "Build to native code only"  ON)
if(DILL_IGNORE_NATIVE)
  set(NATIVE_CG FALSE)
  set(NATIVE_ARCH UNSUPPORTED)
endif()
set(LIBFFI_INTERNAL OFF)
if (DILL_NATIVE_ONLY)
  set(LIBFFI_FOUND FALSE)
else()
  find_package(LibFFI)
endif()
if(LIBFFI_FOUND)
  message(STATUS "Enabling emulation")
  set(EMULATION_POSSIBLE TRUE)
elseif((DILL_IGNORE_NATIVE OR (NATIVE_ARCH STREQUAL "UNSUPPORTED")) AND NOT DILL_NATIVE_ONLY)
  find_program (AUTOCONF autoconf)
  find_program (AUTOMAKE automake)
  if ((AUTOCONF STREQUAL "AUTOCONF-NOTFOUND") OR (AUTOMAKE STREQUAL "AUTOMAKE-NOTFOUND"))
    if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
      # I am top-level project.
      message(FATAL_ERROR "DILL_IGNORE_NATIVE set or native architecture unsupported, but autoconf or automake not found, so unable to build libffi")
    else()
      # I am called from other project with add_subdirectory().
      message(STATUS "DILL_IGNORE_NATIVE set or native architecture unsupported, but autoconf or automake not found, so unable to build libffi")
      return()
    endif()
  endif()
  set(LIBFFI_INTERNAL ON)
  message(STATUS "Using private copy of libffi")
  include(ExternalProject)
  if(BUILD_SHARED_LIBS)
    set(libffi_shared_args --enable-shared --disable-static)
  else()
    set(libffi_shared_args --disable-shared --enable-static)
  endif()
  ExternalProject_Add(libffi
    PREFIX ${CMAKE_BINARY_DIR}/libffi
    INSTALL_DIR ${CMAKE_BINARY_DIR}/libffi/install
    URL https://github.com/libffi/libffi/archive/master.zip
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND
      ./autogen.sh &&
      ./configure
        ${libffi_shared_args}
        --prefix=<INSTALL_DIR>
        --includedir=<INSTALL_DIR>/include
        --libdir=<INSTALL_DIR>/lib
        --disable-multi-os-directory
        --disable-docs
  )
  ExternalProject_Get_property(libffi INSTALL_DIR)
  set(LIBFFI_PREFIX "${INSTALL_DIR}")
  install(DIRECTORY ${LIBFFI_PREFIX}/lib/ DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # This directory has to exist at configure time or we get a buggy CMake error
  file(MAKE_DIRECTORY "${LIBFFI_PREFIX}/include")

  set(LIBFFI_INCLUDE_DIRS ${LIBFFI_PREFIX}/include)
  if(BUILD_SHARED_LIBS)
    set(LIBFFI_LIBRARIES ${LIBFFI_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ffi${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(LIBFFI_LIBRARIES ${LIBFFI_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ffi${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()
  set(PC_LIBFFI_FOUND IGNORE)
  find_package(LibFFI)
  add_dependencies(libffi::libffi libffi)

  if(NOT BUILD_SHARED_LIBS)
    install(FILES ${PROJECT_SOURCE_DIR}/cmake/FindLibFFI.cmake
      DESTINATION "${CMAKE_INSTALL_CMAKEDIR}" COMPONENT ${DILL_HEADER_COMPONENT})
  endif()

  message(STATUS "Enabling emulation")
  set(EMULATION_POSSIBLE TRUE)
else()
  message(STATUS "Disabling emulation")
  set(EMULATION_POSSIBLE FALSE)
endif()

find_package(Perl REQUIRED)

option(DILL_ENABLE_DISASSEMBLY
  "enable binutils-based disassembly (default is OFF)" OFF)

set(ARCHITECTURES sparc ppc64le virtual x86 x86_64 ia64 arm5 arm6 arm8)
foreach(_arch ${ARCHITECTURES})
  add_custom_command(
    OUTPUT dill_${_arch}.c
    COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${_arch}.ops
    DEPENDS ${_arch}.ops
  ) 
  list(APPEND arch_files dill_${_arch}.c ${_arch}.c)
  if(NOT _arch STREQUAL "virtual")
    list(APPEND arch_files ${_arch}_rt.c)
  endif()
endforeach()

add_custom_command(
  OUTPUT dill.h dill.c
  COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/base.ops ${BASE_OPS_ARG}
  DEPENDS base.ops
) 

message (STATUS "base.ops argument is ${BASE_OPS_ARG}")

option(DILL_MULTI_TARGET  "Build all available target architectures" OFF)

if(DILL_MULTI_TARGET)
  message(STATUS "arch files is ${arch_files}")
  list(APPEND SRC_list ${arch_files})
elseif(NOT ${NATIVE_ARCH} STREQUAL "UNSUPPORTED")
  list(APPEND SRC_list dill_${ARCH_FILE}.c ${ARCH_FILE}.c ${ARCH_FILE}_rt.c)
endif()

# The MULTI_TARGET variable is only used to generate the preprocessor symbol
# in config.h.  Use the DIL_MULTI_CONFIG option variable for all associated
# logic in the build.
set(MULTI_TARGET ${DILL_MULTI_TARGET})

if(DILL_MULTI_TARGET OR NATIVE_CG OR EMULATION_POSSIBLE)
  if(EMULATION_POSSIBLE)
    list(APPEND TARGET_DEP_TGT libffi::libffi)
    list(APPEND PKG_DEP_PKG libffi)
  endif()
  if(NOT DILL_MULTI_TARGET AND NOT NATIVE_CG)
    message(STATUS
      "Configure has detected no native dynamic code generation support 
   for this architecture (\"${NATIVE_ARCH}\").  However, the libffi
   library was found, so a library configured for emulation will be built."
    )
    set(EMULATION_ONLY TRUE)
    set(BUILD_EMULATOR TRUE)
  endif()
else()
  if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    # I am top-level project.
    message(FATAL_ERROR
	"Configure has detected no native dynamic code generation support 
	 for this architecture (\"${NATIVE_ARCH}\"), -DDILL_MULTI_TARGET=ON 
	 was not specified, and no emulation is possible (libffi library 
	 not found)  NO LIBRARY WILL BE BUILT" )
  else()
    message(STATUS
	"No native dynamic code generation support for this architecture 
	(\"${CMAKE_SYSTEM_PROCESSOR}\"), -DDILL_MULTI_TARGET=ON was not specified, and no emulation 
	is possible (libffi library not found) NO DILL LIBRARY WILL BE BUILT, but this is not a failure")
    return()
  endif()
endif()

include(CheckIncludeFiles)
check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(stdarg.h STDC_HEADERS)
check_include_files(malloc.h HAVE_MALLOC_H)
check_include_files(memory.h HAVE_MEMORY_H)
check_include_files(sys/mman.h HAVE_SYS_MMAN_H)
include(CheckSymbolExists)
check_symbol_exists(__clear_cache "" CLEAR_CACHE_DEFINED)

set(NO_DISASSEMBLER TRUE)
if(DILL_ENABLE_DISASSEMBLY)
  find_package(BinUtils REQUIRED)
  #list(APPEND TARGET_DEP_TGT binutils::opcodes)
  list(APPEND TARGET_DEP_TGT binutils::bfd)
  list(APPEND PKG_DEP_LIBS ${BinUtils_LIBRARIES})

  set(HAVE_DIS_ASM_H TRUE)
  set(NO_DISASSEMBLER FALSE)

  include(CheckCSourceCompiles)
  include(CMakePushCheckState)
  cmake_push_check_state(RESET)

  set(CMAKE_REQUIRED_INCLUDES ${BinUtils_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${BinUtils_LIBRARIES})
  check_c_source_compiles([=[
#include "dis-asm.h"
void main()
{
  struct disassemble_info info;
  INIT_DISASSEMBLE_INFO(info, stdout, fprintf);
}
]=]
    INIT_DISASSEMBLE_INFO_THREE_ARG
  )

  include(CheckFunctionExists)
  check_function_exists(print_insn_sparc          HAVE_PRINT_INSN_SPARC)
  check_function_exists(print_insn_big_powerpc    HAVE_PRINT_INSN_BIG_POWERPC)
  check_function_exists(print_insn_little_powerpc HAVE_PRINT_INSN_LITTLE_POWERPC)
  check_function_exists(print_insn_ia64           HAVE_PRINT_INSN_IA64)
  check_function_exists(print_insn_i386           HAVE_PRINT_INSN_I386)
  check_function_exists(print_insn_arm            HAVE_PRINT_INSN_ARM)
  check_function_exists(print_insn_little_arm     HAVE_PRINT_INSN_LITTLE_ARM)

  cmake_pop_check_state()
endif()

set(DILL_LIBRARY_PREFIX "" CACHE STRING
  "Prefix to prepend to the output library name")
mark_as_advanced(DILL_LIBRARY_PREFIX)

add_library(dill ${SRC_list})
set_target_properties(dill PROPERTIES
  OUTPUT_NAME ${DILL_LIBRARY_PREFIX}dill
  VERSION ${DILL_VERSION}
  SOVERSION ${DILL_VERSION_MAJOR})
add_library(dill::dill ALIAS dill)

target_include_directories(dill
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${TARGET_DEP_INC}
)
target_link_libraries(dill PRIVATE ${TARGET_DEP_TGT} ${TARGET_DEP_LIBS})

# Setup pkgconfig
option(DILL_INSTALL_PKGCONFIG "Install Dill pkgconfig files" ON)
mark_as_advanced(DILL_INSTALL_PKGCONFIG)
if(DILL_INSTALL_PKGCONFIG)
  set(_pkg_config_private_libs)
  foreach(L ${PKG_DEP_LIBS})
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
  if(_pkg_config_private_libs)
    list(REMOVE_DUPLICATES _pkg_config_private_libs)
  endif()
  string(REPLACE ";" " " _pkg_config_private_libs "${_pkg_config_private_libs}")
  string(REPLACE ";" " " _pkg_config_req_pkg "${PKG_DEP_PKG}")
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/dill.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/dill.pc
    @ONLY
  )
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/dill.pc
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig" COMPONENT ${DILL_HEADER_COMPONENT})
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/dill-config.in
    ${CMAKE_CURRENT_BINARY_DIR}/dill-config
    @ONLY
  )
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/dill-config
    DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT ${DILL_HEADER_COMPONENT})
endif()

option(DILL_INSTALL_HEADERS "Install Dill header files" ON)
mark_as_advanced(DILL_INSTALL_HEADERS)
if(DILL_INSTALL_HEADERS)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dill.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT ${DILL_HEADER_COMPONENT})
endif()

set(namelink_component_args)
if(NOT CMAKE_VERSION VERSION_LESS 3.12)
  set(namelink_component_args NAMELINK_COMPONENT ${DILL_HEADER_COMPONENT})
endif()
install(TARGETS dill
  # IMPORTANT: Add the dill library to the "export-set"
  EXPORT dill-targets
  RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}"          COMPONENT ${DILL_RUNTIME_COMPONENT}
  LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"          COMPONENT ${DILL_LIBRARY_COMPONENT} ${namelink_component_args}
  ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"          COMPONENT ${DILL_ARCHIVE_COMPONENT}
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dill" COMPONENT ${DILL_HEADER_COMPONENT})

if(${CMAKE_C_COMPILER_ID} MATCHES "Intel") 
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -shared-intel")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# Add all targets to the build-tree export set
export(TARGETS dill NAMESPACE dill::
  FILE "${PROJECT_BINARY_DIR}/dill-targets.cmake")
 
configure_file(dill-config-common.cmake.in
  "${PROJECT_BINARY_DIR}/dill-config-common.cmake" @ONLY)

configure_file(dill-config-version.cmake.in
  "${PROJECT_BINARY_DIR}/dill-config-version.cmake" @ONLY)
 
configure_file(dill-config.cmake.in
  "${PROJECT_BINARY_DIR}/dill-config.cmake" @ONLY)

configure_file(dill-config-install.cmake.in
  "${PROJECT_BINARY_DIR}/dill-config-install.cmake" @ONLY)

# Install the dill-config.cmake and dill-config-version.cmake
install(FILES
  "${PROJECT_BINARY_DIR}/dill-config-common.cmake"
  "${PROJECT_BINARY_DIR}/dill-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_CMAKEDIR}" COMPONENT ${DILL_HEADER_COMPONENT})

install(
  FILES "${PROJECT_BINARY_DIR}/dill-config-install.cmake"
  RENAME dill-config.cmake
  DESTINATION "${CMAKE_INSTALL_CMAKEDIR}" COMPONENT ${DILL_HEADER_COMPONENT})
 
# Install the export set for use with the install-tree
install(EXPORT dill-targets NAMESPACE dill::
  DESTINATION "${CMAKE_INSTALL_CMAKEDIR}" COMPONENT ${DILL_HEADER_COMPONENT})

# Install extra find module dependencies
install(DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR} COMPONENT ${DILL_HEADER_COMPONENT}
  FILES_MATCHING PATTERN "Find*.cmake" PATTERN "CMake*.cmake"
)

# display status message for important variables
include(CTest)
mark_as_advanced(BUILD_TESTING)

if(BUILD_TESTING)
  enable_testing()
  set(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  mark_as_advanced(BUILDNAME)
  if(NATIVE_CG)
    add_subdirectory(tests)
  endif()
  add_subdirectory(vtests)
endif()

configure_file(
  ${PROJECT_SOURCE_DIR}/CTestCustom.cmake.in
  ${PROJECT_BINARY_DIR}/CTestCustom.cmake
  @ONLY
)

option(DILL_QUIET "Suppress summary output" OFF)
if(NOT DILL_QUIET)
  message(STATUS)
  message(STATUS "-----------------------------------------------------------------------------")
  message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
  message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
  message(STATUS "DILL_ENABLE_DISASSEMBLY = ${DILL_ENABLE_DISASSEMBLY}" )
  message(STATUS "DILL_MULTI_TARGET = ${DILL_MULTI_TARGET}" )
  message(STATUS "DILL_IGNORE_NATIVE = ${DILL_IGNORE_NATIVE}" )
  message(STATUS "DILL_NATIVE_ONLY = ${DILL_NATIVE_ONLY}" )
  message(STATUS "BUILD_TESTING = ${BUILD_TESTING}" )
  message(STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
  message(STATUS "Change a value with: cmake -D<Variable>=<Value>" )
  message(STATUS "-----------------------------------------------------------------------------")
endif()
