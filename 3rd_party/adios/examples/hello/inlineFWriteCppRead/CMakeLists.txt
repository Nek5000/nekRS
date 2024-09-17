#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2HelloFWriteCppReadExample)

if(NOT TARGET adios2_core)
  set(_components CXX)

  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  endif()
  if(CMAKE_Fortran_COMPILER_LOADED)
    list(APPEND _components Fortran)
  endif()

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

if(ADIOS2_HAVE_Fortran)
  include(FortranCInterface)
  FortranCInterface_HEADER(InlineExampleFC.h MACRO_NAMESPACE "FC_")
  FortranCInterface_VERIFY(CXX QUIET)

  add_library(adios2_hello_inlineMixedLang inlineMixedLang.cpp)
  target_link_libraries(adios2_hello_inlineMixedLang adios2::cxx11)
  target_include_directories(adios2_hello_inlineMixedLang PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  install(TARGETS adios2_hello_inlineMixedLang RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_hello_inlineMixedLang_f inlineMixedLang.f90)
  target_link_libraries(adios2_hello_inlineMixedLang_f adios2_hello_inlineMixedLang adios2::fortran)
  set_target_properties(adios2_hello_inlineMixedLang_f PROPERTIES LINKER_LANGUAGE Fortran)
  install(TARGETS adios2_hello_inlineMixedLang_f RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
