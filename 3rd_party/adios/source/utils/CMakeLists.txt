#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/verinfo.h.in
  ${PROJECT_BINARY_DIR}/verinfo.h
  @ONLY
)

# BPLS
add_executable(bpls ./bpls/bpls.cpp)
target_link_libraries(bpls
                      PUBLIC adios2_core adios2sys
                      PRIVATE adios2::thirdparty::pugixml $<$<PLATFORM_ID:Windows>:shlwapi>)

target_include_directories(bpls PRIVATE
  ${PROJECT_BINARY_DIR}
  ${PROJECT_SOURCE_DIR}/bindings/C
)

set_property(TARGET bpls PROPERTY OUTPUT_NAME bpls${ADIOS2_EXECUTABLE_SUFFIX})
install(TARGETS bpls EXPORT adios2
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_tools-runtime
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/bpls/bpls.cmake.gen.in
  ${PROJECT_BINARY_DIR}/bpls.cmake.gen
  @ONLY
)
file(GENERATE
  OUTPUT ${PROJECT_BINARY_DIR}/$<CONFIG>/bpls.cmake
  INPUT ${PROJECT_BINARY_DIR}/bpls.cmake.gen
)


# BPSPLIT
#add_executable(bpsplit ./bpsplit/bpsplit.cpp)
#target_link_libraries(bpsplit adios2::cxx11 adios2sys_interface adios2::thirdparty::pugixml)
#install(TARGETS bpsplit EXPORT adios2
#  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#)

# ADIOS_REORGANIZE
set(adios_reorganize_srcs
  adios_reorganize/main.cpp
  adios_reorganize/Reorganize.cpp
  Utils.cpp
  )
if(ADIOS2_HAVE_MGARD)
  set(maybe_mgard mgard::mgard)
else()
  set(maybe_mgard)
endif()
add_executable(adios_reorganize ${adios_reorganize_srcs})
target_link_libraries(adios_reorganize PRIVATE adios2_core ${maybe_mgard})
set_property(TARGET adios_reorganize PROPERTY OUTPUT_NAME adios2_reorganize${ADIOS2_EXECUTABLE_SUFFIX})

if(ADIOS2_HAVE_MPI)
  add_executable(adios_reorganize_mpi ${adios_reorganize_srcs})
  target_link_libraries(adios_reorganize_mpi PRIVATE adios2_core_mpi ${maybe_mgard})
  set_property(TARGET adios_reorganize_mpi PROPERTY OUTPUT_NAME adios2_reorganize_mpi${ADIOS2_EXECUTABLE_SUFFIX})
  set(maybe_adios_reorganize_mpi adios_reorganize_mpi)
else()
  set(maybe_adios_reorganize_mpi)
endif()

install(TARGETS adios_reorganize
  ${maybe_adios_reorganize_mpi}
  EXPORT adios2
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_tools-runtime
)

if(ADIOS2_HAVE_MPI)
  add_subdirectory(adios_iotest)
endif()

if(Python_Interpreter_FOUND)
  add_subdirectory(bp4dbg)
  add_subdirectory(bp5dbg)
  add_subdirectory(adios_campaign_manager)
  install(PROGRAMS adios2_json_pp.py
    RENAME adios2_json_pp
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT adios2_scripts-runtime)
endif()

install(PROGRAMS adios2_deactivate_bp
  RENAME adios2_deactivate_bp
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT adios2_scripts-runtime)

# Simplified wrappers for adios2_reorganize
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/adios_reorganize/adios2_reorganize_wrapper
  ${PROJECT_BINARY_DIR}/adios2_reorganize_wrapper
  @ONLY
)

if(ADIOS2_HAVE_HDF5)
  install(PROGRAMS ${PROJECT_BINARY_DIR}/adios2_reorganize_wrapper
    RENAME bp2h5${ADIOS2_EXECUTABLE_SUFFIX}
    DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
  if(ADIOS2_HAVE_MPI)
    install(PROGRAMS ${PROJECT_BINARY_DIR}/adios2_reorganize_wrapper
      RENAME bp2h5_mpi${ADIOS2_EXECUTABLE_SUFFIX}
      DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
  endif()
endif()
