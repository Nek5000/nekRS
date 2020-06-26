include(ExternalProject)
include(FetchContent)
#set(FETCHCONTENT_QUIET OFF)

# ---------------------------------------------------------
# Download sources
# ---------------------------------------------------------

# Nek5000
# =======

FetchContent_Declare(
  nek5000_content
  GIT_REPOSITORY https://github.com/nek5000/nek5000.git
  GIT_TAG master
  )
FetchContent_GetProperties(nek5000_content)
if (NOT nek5000_content_POPULATED)
  FetchContent_Populate(nek5000_content)
endif()

set(NEK5000_DIR ${nek5000_content_SOURCE_DIR})

# blasLapack
# ==========

set(BLASLAPACK_DIR ${NEK5000_DIR}/3rd_party/blasLapack)

# gslib
# =====

set(GS_DIR ${NEK5000_DIR}/3rd_party/gslib)

FetchContent_Declare(
  gs_content
  URL http://github.com/gslib/gslib/archive/v1.0.5.tar.gz
  URL_HASH SHA1=1b5b28de5b997c3b0893a3b4dcf5cee8614b9f27
  SOURCE_DIR ${GS_DIR}
)
FetchContent_GetProperties(gs_content)
if (NOT gs_content_POPULATED)
  file(RENAME ${GS_DIR}/install temp_gslib_install)
  FetchContent_Populate(gs_content)
  file(RENAME temp_gslib_install ${GS_DIR}/install)
  file(MAKE_DIRECTORY ${GS_DIR}/include)
endif()

# parRSB
# ======

if (PARRSB MATCHES ${NEK5000_PPLIST})
  FetchContent_Declare(
    parrsb_content
    URL https://github.com/Nek5000/parRSB/archive/v0.2.tar.gz
    SOURCE_DIR ${NEK5000_DIR}/3rd_party/parRSB
  )
  FetchContent_GetProperties(parrsb_content)
  if (NOT parrsb_content_POPULATED)
    FetchContent_Populate(parrsb_content)
  endif()
endif()

# ---------------------------------------------------------
# Project configuration
# ---------------------------------------------------------

# ExternalProject_Add has some peculiar rules about quoting arguments,
# so we use the helper script (run_config.sh) to set the environment
# variables based on the command-line arguments
ExternalProject_Add(
  nek5000_project
  SOURCE_DIR ${NEK5000_DIR}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND 
    ${CMAKE_CURRENT_LIST_DIR}/run_nekconfig.sh 
    "CC=${CMAKE_C_COMPILER}" 
    "CFLAGS=${CMAKE_C_FLAGS}" 
    "FC=${CMAKE_Fortran_COMPILER}"
    "FFLAGS=${CMAKE_Fortran_FLAGS}"
    "NEK5000_DIR=${NEK5000_DIR}"
  INSTALL_COMMAND ""
  USES_TERMINAL_BUILD on
)

add_library(gs STATIC IMPORTED)
target_include_directories(gs INTERFACE ${GS_DIR}/src ${GS_DIR}/include)
set_target_properties(gs PROPERTIES IMPORTED_LOCATION ${GS_DIR}/lib/libgs.a)
add_dependencies(gs nek5000_project)

add_library(blasLapack STATIC IMPORTED)
set_target_properties(blasLapack PROPERTIES IMPORTED_LOCATION ${BLASLAPACK_DIR}/libblasLapack.a)
add_dependencies(blasLapack nek5000_project)
