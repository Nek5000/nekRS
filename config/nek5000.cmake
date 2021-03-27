set(NEK5000_GS_VERSION "1.0.7")
set(PARRSB_VERSION "0.7")

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

if (USE_PARRSB)
  install(FILES ${NEK5000_SOURCE_DIR}/core/PARALLEL.dprocmap DESTINATION ${NEK5000_SOURCE_DIR}/core RENAME "PARALLEL")
else()
  install(FILES ${NEK5000_SOURCE_DIR}/core/PARALLEL.default DESTINATION ${NEK5000_SOURCE_DIR}/core RENAME "PARALLEL")
endif()

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
#set(NEK5000_GS_SOURCE_DIR ${NEK5000_GS_DIR}/gslib/src)
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

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(FPIC_FLAG "-fPIC")

  set(MCMODEL_FLAG "-mcmodel=medium")
  CHECK_C_COMPILER_FLAG("${MCMODEL_FLAG}" COMPILER_C_SUPPORTS_MCMODEL_MEDIUM)
  if(NOT COMPILER_C_SUPPORTS_MCMODEL_MEDIUM) 
    set(MCMODEL_FLAG "")
  endif()

  CHECK_C_COMPILER_FLAG("${MCMODEL_FLAG} ${MCMODEL_FLAG}" COMPILER_C_SUPPORTS_MCMODEL_MEDIUM_FPIC)
  if(NOT COMPILER_C_SUPPORTS_MCMODEL_MEDIUM_FPIC) 
    set(MCMODEL_FLAG "")
  endif()
endif()

include(CheckFortranCompilerFlag)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  CHECK_Fortran_COMPILER_FLAG("-fcray-pointer" COMPILER_Fortran_SUPPORTS_CRAYPTR)
  if(COMPILER_Fortran_SUPPORTS_CRAYPTR)
    set(CRAYPTR_FLAG "-fcray-pointer")
  else()
    MESSAGE(FATAL_ERROR "Compiler does not support Cray pointers!")
  endif()
endif()

ExternalProject_Add(
  nek5000_deps
  SOURCE_DIR ${NEK5000_SOURCE_DIR}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND 
    ${CMAKE_CURRENT_LIST_DIR}/run_nekconfig.sh 
    "CC=${CMAKE_C_COMPILER}" 
    "CFLAGS=${EXTERNAL_C_FLAGS} ${FPIC_FLAG} ${MCMODEL_FLAG}"
    "FC=${CMAKE_Fortran_COMPILER}"
    "FFLAGS=${EXTERNAL_Fortran_FLAGS} ${FPIC_FLAG} ${MCMODEL_FLAG} ${CRAYPTR_FLAG} -I../../"
    "NEK5000_SOURCE_DIR=${NEK5000_SOURCE_DIR}"
    "PPLIST=${NEK5000_PPLIST}"
  INSTALL_COMMAND ""
  USES_TERMINAL_BUILD on
)

unset(MCMODEL_FLAG)
unset(CRAYPTR_FLAG)
unset(FPIC_FLAG)

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

