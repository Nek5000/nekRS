include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

set(OCCA_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/occa)

# If OCCA is not in source tree, download it in build directory
if (EXISTS ${OCCA_SOURCE_DIR}) 
  message(STATUS "Using OCCA source in ${OCCA_SOURCE_DIR}")
  add_subdirectory(${OCCA_SOURCE_DIR})
else()
  FetchContent_Declare(
    occa_content
    GIT_REPOSITORY https://github.com/libocca/occa.git
    GIT_TAG master)
  FetchContent_MakeAvailable(occa_content)
  FetchContent_GetProperties(occa_content)
  set(OCCA_SOURCE_DIR ${occa_content_SOURCE_DIR})
endif()

#set(occa_content_BUILD_TESTS OFF)

