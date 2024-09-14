if(NOT FFS_LIBRARIES_ONLY)
  add_subdirectory(progs)
endif()

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
