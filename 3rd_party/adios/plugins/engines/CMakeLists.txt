#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if(ADIOS2_HAVE_Catalyst)
  add_library(ParaViewADIOSInSituEngine MODULE
    ParaViewFidesEngine.cpp)
  target_link_libraries(ParaViewADIOSInSituEngine
    PUBLIC
      adios2_core
    PRIVATE
      catalyst::catalyst)
  install(TARGETS ParaViewADIOSInSituEngine  EXPORT adios2Exports
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_core-runtime
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT adios2_core-libraries NAMELINK_COMPONENT adios2_core-development
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT adios2_core-development
  )
endif()
