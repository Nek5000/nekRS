set(GS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/gslib/)
set(GS_LIB 3rd_party/gslib/libgs.a)
set(OGS_SOURCE_DIR ${GS_SOURCE_DIR}/ogs)

# =============================================================
# Build gslib
# =============================================================

# Copy source into CMake build dir.  We do this since gslib is built in-source,
# and we want to keep source tree clean.
FetchContent_Declare(
    gs_content
    URL ${GS_SOURCE_DIR} 
)
FetchContent_GetProperties(gs_content)
if (NOT gs_content_POPULATED)
    FetchContent_Populate(gs_content)
endif()
set(GS_SOURCE_DIR ${gs_content_SOURCE_DIR})

# Build gslib
ExternalProject_Add(
        gs_build
        SOURCE_DIR ${GS_SOURCE_DIR}
        BUILD_IN_SOURCE on
        CONFIGURE_COMMAND ""
        BUILD_COMMAND 
          ${CMAKE_CURRENT_LIST_DIR}/run_gs_install.sh
          "CC=${CMAKE_C_COMPILER}"
          "CFLAGS=-fPIC ${EXTERNAL_C_FLAGS}"
        INSTALL_COMMAND ""
        USES_TERMINAL_BUILD on
)

# Target for libraries
add_library(gs STATIC IMPORTED)
set_target_properties(gs PROPERTIES IMPORTED_LOCATION ${GS_SOURCE_DIR}/lib/libgs.a)
target_include_directories(gs INTERFACE ${GS_SOURCE_DIR}/src)
add_dependencies(gs gs_build)

# =============================================================
# OGS Sources
# =============================================================

set(OGS_SOURCES
        ${OGS_SOURCE_DIR}/src/ogsGather.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherMany.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherScatter.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherScatterMany.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherScatterVec.cpp
        ${OGS_SOURCE_DIR}/src/ogsGatherVec.cpp
        ${OGS_SOURCE_DIR}/src/ogsHostGather.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherMany.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherScatter.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherScatterMany.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherScatterVec.c
        ${OGS_SOURCE_DIR}/src/ogsHostGatherVec.c
        ${OGS_SOURCE_DIR}/src/ogsHostScatter.c
        ${OGS_SOURCE_DIR}/src/ogsHostScatterMany.c
        ${OGS_SOURCE_DIR}/src/ogsHostScatterVec.c
        ${OGS_SOURCE_DIR}/src/ogsHostSetup.c
        ${OGS_SOURCE_DIR}/src/ogsKernels.cpp
        ${OGS_SOURCE_DIR}/src/ogsMappedAlloc.cpp
        ${OGS_SOURCE_DIR}/src/ogsScatter.cpp
        ${OGS_SOURCE_DIR}/src/ogsScatterMany.cpp
        ${OGS_SOURCE_DIR}/src/ogsScatterVec.cpp
        ${OGS_SOURCE_DIR}/src/ogsSetup.cpp
        ${OGS_SOURCE_DIR}/src/oogs.cpp)

set(file_pattern "\.cu$|\.hip$|\.okl$|\.c$|\.hpp$|\.tpp$|\.h$$")

install(DIRECTORY
        ${GS_SOURCE_DIR}/src/
        DESTINATION gslib
        FILES_MATCHING REGEX "\.h$")
install(DIRECTORY
        ${OGS_SOURCE_DIR}/include
        ${OGS_SOURCE_DIR}/okl
        DESTINATION gatherScatter
        FILES_MATCHING REGEX ${file_pattern})
install(FILES ${OGS_SOURCE_DIR}/ogs.hpp DESTINATION gatherScatter)
