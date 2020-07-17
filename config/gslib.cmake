set(LIBP_GS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/gslib/)
set(LIBP_GS_LIB 3rd_party/gslib/libgs.a)
set(OGS_SOURCE_DIR ${LIBP_GS_SOURCE_DIR}/ogs)

# =============================================================
# Build gslib
# =============================================================

# Copy source into CMake build dir.  We do this since gslib is built in-source,
# and we want to keep source tree clean.
FetchContent_Declare(
    libp_gs_content
    URL ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/gslib
)
FetchContent_GetProperties(libp_gs_content)
if (NOT libp_gs_content)
    FetchContent_Populate(libp_gs_content)
endif()
set(LIBP_GS_SOURCE_DIR ${libp_gs_content_SOURCE_DIR})

# Build gslib
ExternalProject_Add(
        libp_gs_build
        SOURCE_DIR ${LIBP_GS_SOURCE_DIR}
        BUILD_IN_SOURCE on
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${CMAKE_CURRENT_LIST_DIR}/run_libp_gs_install.sh "CC=${CMAKE_C_COMPILER}"
        INSTALL_COMMAND ""
        USES_TERMINAL_BUILD on
)

# Target for libraries
add_library(libp_gs STATIC IMPORTED)
set_target_properties(libp_gs PROPERTIES IMPORTED_LOCATION ${LIBP_GS_SOURCE_DIR}/lib/libgs.a)
target_include_directories(libp_gs INTERFACE ${LIBP_GS_SOURCE_DIR}/src)
add_dependencies(libp_gs libp_gs_build)

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

install(DIRECTORY
        ${OGS_SOURCE_DIR}/include
        ${OGS_SOURCE_DIR}/okl
        DESTINATION gatherScatter
        FILES_MATCHING REGEX ${file_pattern})
install(FILES ${OGS_SOURCE_DIR}/ogs.hpp DESTINATION gatherScatter)




