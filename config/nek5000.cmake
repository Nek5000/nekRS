message(CHECK_START "Checking for a supported Fortran compiler")
string(COMPARE EQUAL "${CMAKE_Fortran_COMPILER_ID}" "GNU" USING_GNU)
string(COMPARE EQUAL "${CMAKE_Fortran_COMPILER_ID}" "IntelLLVM" USING_INTEL_LLVM)
string(COMPARE EQUAL "${CMAKE_Fortran_COMPILER_ID}" "NVHPC" USING_NVHPC)

if(USING_GNU OR USING_INTEL_LLVM OR USING_NVHPC)
  message(CHECK_PASS "Found the ${CMAKE_Fortran_COMPILER_ID} Fortran compiler")
else()
  message(FATAL_ERROR "No supported Fortran compiles found to build Nek5000 interface!")
endif()

if (${NEK5000_PPLIST} MATCHES "PARRSB")
  set(USE_PARRSB on)
endif()

# ---------------------------------------------------------
# Download sources
# ---------------------------------------------------------

# Nek5000
# =======

  # Since Nek5000 is compiled in-source, we copy it to the build 
FetchContent_Declare(
  nek5000_content
  URL ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/nek5000)
FetchContent_GetProperties(nek5000_content)
if (NOT nek5000_content_POPULATED)
  FetchContent_Populate(nek5000_content)
endif()

set(NEK5000_SOURCE_DIR ${nek5000_content_SOURCE_DIR})

# blasLapack
# ==========

set(BLASLAPACK_DIR ${NEK5000_SOURCE_DIR}/3rd_party/blasLapack)

# gslib
# =====

set(NEK5000_GS_SUBTREE ${CMAKE_CURRENT_LIST_DIR}/../3rd_party/nek5000_gslib)
set(NEK5000_GS_DIR ${NEK5000_SOURCE_DIR}/3rd_party/gslib/gslib)

FetchContent_Declare(
  nek5000_gs_content
  URL ${NEK5000_GS_SUBTREE}
  SOURCE_DIR ${NEK5000_GS_DIR}
)
if (NOT nek5000_gs_content_POPULATED)
  FetchContent_Populate(nek5000_gs_content)
endif()

# ./build/_deps/nek5000_content-src/3rd_party/gslib/lib/libgs.a
set(NEK5000_GS_INCLUDE_DIR ${NEK5000_GS_DIR}/../include)
set(NEK5000_GS_LIB_DIR ${NEK5000_GS_DIR}/../lib)


# parRSB
# ======

set(PARRSB_SUBTREE ${CMAKE_CURRENT_LIST_DIR}/../3rd_party/nek5000_parRSB)
set(PARRSB_DIR ${NEK5000_SOURCE_DIR}/3rd_party/parRSB/parRSB)

FetchContent_Declare(
  parrsb_content
  URL ${PARRSB_SUBTREE}
  SOURCE_DIR ${PARRSB_DIR}
)
FetchContent_GetProperties(parrsb_content)
if (NOT parrsb_content_POPULATED)
  # Moves `install` script so it doesn't get clobbered when 
  # populating content.  
  #file(RENAME ${PARRSB_DIR}/install temp_parrsb_install)
  FetchContent_Populate(parrsb_content)
  #file(RENAME temp_parrsb_install ${PARRSB_DIR}/install)
endif()

set(PARRSB_INCLUDE_DIR ${PARRSB_DIR}/../include)
set(PARRSB_LIB_DIR ${PARRSB_DIR}/../lib)

# ---------------------------------------------------------
# Build Nek5000 dependencies
# ---------------------------------------------------------

# ExternalProject_Add has some peculiar rules about quoting arguments,
# so we use the helper script (run_config.sh) to set the environment
# variables based on the command-line arguments

set(FPIC_FLAG "-fPIC")

set(PPC FALSE)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(powerpc|ppc)64")
set(PPC TRUE)
endif()

set(MCMODEL_FLAG "-mcmodel=medium")
if(PPC)
  set(MCMODEL_FLAG "-mcmodel=large")
elseif(APPLE)
  set(MCMODEL_FLAG "")
endif()

include(CheckCCompilerFlag)
if(NOT MCMODEL_FLAG STREQUAL "")
  CHECK_C_COMPILER_FLAG("${MCMODEL_FLAG}" COMPILER_C_SUPPORTS_MCMODEL_MEDIUM)
  if(NOT COMPILER_C_SUPPORTS_MCMODEL_MEDIUM AND NOT PPC)
    set(MCMODEL_FLAG "-mcmodel=large")
    CHECK_C_COMPILER_FLAG("${MCMODEL_FLAG}" COMPILER_C_SUPPORTS_MCMODEL_LARGE)
    if(NOT COMPILER_C_SUPPORTS_MCMODEL_LARGE) 
      set(MCMODEL_FLAG "")
    endif()
  endif()
endif()

if(NOT MCMODEL_FLAG STREQUAL "")
  CHECK_C_COMPILER_FLAG("${FPIC_FLAG} ${MCMODEL_FLAG}" COMPILER_C_SUPPORTS_MCMODEL_FPIC)
  if(NOT COMPILER_C_SUPPORTS_MCMODEL_FPIC)
   set(MCMODEL_FLAG "")
   endif()
endif()

include(CheckFortranCompilerFlag)
if(MCMODEL_FLAG STREQUAL "-mcmodel=medium" AND USING_GNU)
  set(LARGE_DATA_THRES_FLAG "-mlarge-data-threshold=0")
  CHECK_Fortran_COMPILER_FLAG("${MCMODEL_FLAG} ${LARGE_DATA_THRES_FLAG}" COMPILER_Fortran_SUPPORTS_LARGE_DATA_THRES)
  if(NOT COMPILER_Fortran_SUPPORTS_LARGE_DATA_THRES) 
    set(LARGE_DATA_THRES_FLAG "")
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CRAYPTR_FLAG "-fcray-pointer")
endif()

string(REGEX REPLACE "-O[2,3,fast]" "" _EXTERNAL_C_FLAGS ${EXTERNAL_C_FLAGS})
string(REGEX REPLACE "-O[2,3,fast]" "" _EXTERNAL_Fortran_FLAGS ${EXTERNAL_Fortran_FLAGS})

ExternalProject_Add(
  nek5000_deps
  SOURCE_DIR ${NEK5000_SOURCE_DIR}
  CONFIGURE_COMMAND
    ${CMAKE_CURRENT_LIST_DIR}/run_nekconfig.sh 
    "LDFLAGS=${BSYMBOLIC_FLAG}" 
    "CC=${CMAKE_C_COMPILER}" 
    "CFLAGS=${_EXTERNAL_C_FLAGS} ${FPIC_FLAG} ${MCMODEL_FLAG} ${LARGE_DATA_THRES_FLAG}"
    "FC=${CMAKE_Fortran_COMPILER}"
    "FFLAGS=${_EXTERNAL_Fortran_FLAGS} ${FPIC_FLAG} ${MCMODEL_FLAG} ${LARGE_DATA_THRES_FLAG} ${CRAYPTR_FLAG}"
    "NEK5000_SOURCE_DIR=${NEK5000_SOURCE_DIR}"
    "PPLIST=${NEK5000_PPLIST}"

  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  USES_TERMINAL_BUILD on
)

unset(MCMODEL_FLAG)
unset(CRAYPTR_FLAG)
unset(FPIC_FLAG)
unset(PPC)

add_library(nek5000_gs STATIC IMPORTED)
set_target_properties(nek5000_gs PROPERTIES IMPORTED_LOCATION ${NEK5000_GS_LIB_DIR}/libgs.a)
target_include_directories(nek5000_gs INTERFACE ${NEK5000_GS_SOURCE_DIR} ${NEK5000_GS_INCLUDE_DIR})

add_dependencies(nek5000_gs nek5000_deps)

add_library(blasLapack STATIC IMPORTED)
set_target_properties(blasLapack PROPERTIES IMPORTED_LOCATION ${BLASLAPACK_DIR}/libblasLapack.a)
add_dependencies(blasLapack nek5000_deps)

if (${USE_PARRSB})
  add_library(parRSB STATIC IMPORTED)
  set_target_properties(parRSB PROPERTIES IMPORTED_LOCATION ${PARRSB_DIR}/lib/libparRSB.a)
  add_dependencies(parRSB nek5000_deps)
endif()

# ---------------------------------------------------------
# Install
# ---------------------------------------------------------


install(DIRECTORY ${NEK5000_SOURCE_DIR}/core DESTINATION nek5000
  PATTERN "*"
  PATTERN "mkuserfile" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(PROGRAMS ${NEK5000_SOURCE_DIR}/bin/nekconfig DESTINATION nek5000/bin)

ExternalProject_Get_property(nek5000_deps BINARY_DIR)
install(FILES ${BINARY_DIR}/makefile.template DESTINATION nek5000)

install(FILES ${NEK5000_GS_LIB_DIR}/libgs.a DESTINATION nek5000/3rd_party/gslib/lib)
install(DIRECTORY ${NEK5000_GS_INCLUDE_DIR} DESTINATION nek5000/3rd_party/gslib)

install(FILES ${BLASLAPACK_DIR}/libblasLapack.a DESTINATION nek5000/3rd_party/blasLapack)

if (${USE_PARRSB})
  install(FILES ${PARRSB_LIB_DIR}/libparRSB.a DESTINATION nek5000/3rd_party/parRSB/lib)
  install(DIRECTORY ${PARRSB_INCLUDE_DIR} DESTINATION nek5000/3rd_party/parRSB)
endif()
