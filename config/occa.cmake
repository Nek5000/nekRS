include(FetchContent)
#set(FETCHCONTENT_QUIET OFF)

FetchContent_Declare(
  occa_content
  GIT_REPOSITORY https://github.com/libocca/occa.git
  GIT_TAG master)

#set(occa_content_BUILD_TESTS OFF)
FetchContent_MakeAvailable(occa_content)

