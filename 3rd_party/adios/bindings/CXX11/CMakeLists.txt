#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

add_library(adios2_cxx11
  adios2/cxx11/ADIOS.cpp
  adios2/cxx11/ADIOSView.h
  adios2/cxx11/Attribute.cpp
  adios2/cxx11/Engine.cpp
  adios2/cxx11/Engine.tcc
  adios2/cxx11/IO.cpp
  adios2/cxx11/IO.tcc
  adios2/cxx11/Operator.cpp
  adios2/cxx11/Query.cpp
  adios2/cxx11/Types.cpp
  adios2/cxx11/Types.tcc
  adios2/cxx11/Variable.cpp
  adios2/cxx11/VariableNT.cpp
  adios2/cxx11/fstream/ADIOS2fstream.cpp
  adios2/cxx11/Group.cpp
  adios2/cxx11/Group.tcc
)

set_property(TARGET adios2_cxx11 PROPERTY EXPORT_NAME cxx11)
set_property(TARGET adios2_cxx11 PROPERTY OUTPUT_NAME adios2${ADIOS2_LIBRARY_SUFFIX}_cxx11)
target_link_libraries(adios2_cxx11 PRIVATE adios2_core)
target_compile_features(adios2_cxx11 INTERFACE ${ADIOS2_CXX11_FEATURES})

target_include_directories(adios2_cxx11
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  INTERFACE
    $<BUILD_INTERFACE:$<TARGET_PROPERTY:adios2_core,INTERFACE_INCLUDE_DIRECTORIES>>
)

add_library(adios2::cxx11 ALIAS adios2_cxx11)

if (ADIOS2_HAVE_Derived_Variable)
    target_sources(adios2_cxx11 PRIVATE
      adios2/cxx11/VariableDerived.cpp
  )
endif()

if(ADIOS2_HAVE_MPI)
  add_library(adios2_cxx11_mpi
    adios2/cxx11/ADIOSMPI.cpp
    adios2/cxx11/IOMPI.cpp
    adios2/cxx11/fstream/ADIOS2fstreamMPI.cpp
  )
  set_property(TARGET adios2_cxx11_mpi PROPERTY EXPORT_NAME cxx11_mpi)
  set_property(TARGET adios2_cxx11_mpi PROPERTY OUTPUT_NAME adios2${ADIOS2_LIBRARY_SUFFIX}_cxx11_mpi)
  target_link_libraries(adios2_cxx11_mpi PUBLIC adios2_cxx11 PRIVATE adios2_core_mpi PUBLIC MPI::MPI_C)
  set(maybe_adios2_cxx11_mpi adios2_cxx11_mpi)
  target_compile_definitions(adios2_cxx11_mpi INTERFACE ADIOS2_USE_MPI)
  add_library(adios2::cxx11_mpi ALIAS adios2_cxx11_mpi)
else()
  set(maybe_adios2_cxx11_mpi)
endif()

# Set library version information
set_target_properties(
  adios2_cxx11 ${maybe_adios2_cxx11_mpi}
  PROPERTIES
  VERSION ${ADIOS2_LIBRARY_VERSION}
  SOVERSION ${ADIOS2_LIBRARY_SOVERSION}
)

install(TARGETS adios2_cxx11 ${maybe_adios2_cxx11_mpi} EXPORT adios2CXX11Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_cxx11-runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT adios2_cxx11-libraries NAMELINK_COMPONENT adios2_cxx11-development
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT adios2_cxx11-development
)

install(
  FILES adios2.h 
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT adios2_cxx11-development
)

install(
  FILES adios2/cxx11/fstream/ADIOS2fstream.h 
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/adios2/cxx11/fstream
  COMPONENT adios2_cxx11-development
)

if (ADIOS2_HAVE_Derived_Variable)
    install(
        FILES adios2/cxx11/VariableDerived.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/adios2/cxx11
        COMPONENT adios2_cxx11-development
    )
endif()

install(
  FILES adios2/cxx11/ADIOS.h
        adios2/cxx11/ADIOS.inl
        adios2/cxx11/ADIOSView.h
        adios2/cxx11/IO.h
        adios2/cxx11/Group.h
        adios2/cxx11/Variable.h
        adios2/cxx11/VariableNT.h
        adios2/cxx11/Attribute.h
        adios2/cxx11/Engine.h
        adios2/cxx11/KokkosView.h
        adios2/cxx11/Operator.h
        adios2/cxx11/Query.h
        adios2/cxx11/Types.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/adios2/cxx11
  COMPONENT adios2_cxx11-development
)

export(EXPORT adios2CXX11Exports
  FILE ${ADIOS2_BINARY_DIR}/adios2-cxx11-targets.cmake
  NAMESPACE adios2::
)
install(EXPORT adios2CXX11Exports
  FILE adios2-cxx11-targets.cmake
  NAMESPACE adios2::
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT adios2_cxx11-development
)
