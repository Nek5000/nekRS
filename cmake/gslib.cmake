set(GS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/gslib/)
set(GS_LIB 3rd_party/gslib/libgs.a)

# Copy source into CMake build dir.  We do this since gslib is built in-source,
# and we want to keep source tree clean.
FetchContent_Declare(
    gs_content
    URL ${GS_SOURCE_DIR} 
)
FetchContent_GetProperties(gs_content)
if (NOT gs_content_POPULATED)
    FetchContent_MakeAvailable(gs_content)
endif()
set(GS_SOURCE_DIR ${gs_content_SOURCE_DIR})

# Build gslib
ExternalProject_Add(
        gs_build
        SOURCE_DIR ${GS_SOURCE_DIR}
        BUILD_IN_SOURCE on
        CONFIGURE_COMMAND ""
        BUILD_COMMAND "" 
        INSTALL_COMMAND cd ${GS_SOURCE_DIR} && $(MAKE) CC=${CMAKE_C_COMPILER} "CFLAGS=-fPIC ${EXTERNAL_C_FLAGS}" install
        USES_TERMINAL_BUILD on
)

# Target for libraries
add_library(gs STATIC IMPORTED)
add_dependencies(gs gs_build)
set_target_properties(gs PROPERTIES IMPORTED_LOCATION ${GS_SOURCE_DIR}/build/lib/libgs.a)
file(MAKE_DIRECTORY ${GS_SOURCE_DIR}/build/include)
target_include_directories(gs INTERFACE ${GS_SOURCE_DIR}/build/include)

set(file_pattern "\.cu$|\.hip$|\.okl$|\.c$|\.hpp$|\.tpp$|\.h$$")

install(DIRECTORY
        ${GS_SOURCE_DIR}/build/include 
        DESTINATION gslib
        FILES_MATCHING REGEX "\.h$")
