#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if (NOT ADIOS2_USE_PIP)
  add_executable(adios2_remote_server ./remote_server.cpp remote_common.cpp)

  target_link_libraries(adios2_remote_server
                        PUBLIC EVPath::EVPath adios2_core adios2sys
                        PRIVATE adios2::thirdparty::pugixml $<$<PLATFORM_ID:Windows>:shlwapi>)

  target_include_directories(adios2_remote_server PRIVATE ${PROJECT_BINARY_DIR})

  set_property(TARGET adios2_remote_server PROPERTY OUTPUT_NAME adios2_remote_server${ADIOS2_EXECUTABLE_SUFFIX})
  install(TARGETS adios2_remote_server EXPORT adios2
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_tools-runtime
  )
endif ()
