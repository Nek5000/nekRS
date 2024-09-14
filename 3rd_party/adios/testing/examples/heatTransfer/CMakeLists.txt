#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

include(${CMAKE_CURRENT_SOURCE_DIR}/TestBPFileMx1.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/TestBPFileMxM.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/TestBPFileMxN.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/TestInlineMxM.cmake)

if(ADIOS2_HAVE_ZFP)
  #include(${CMAKE_CURRENT_SOURCE_DIR}/TestBPFileMxM_zfp.cmake)
  #include(${CMAKE_CURRENT_SOURCE_DIR}/TestBPFileMxN_zfp.cmake)
  include(${CMAKE_CURRENT_SOURCE_DIR}/TestBPFileMx1_zfp.cmake)
endif()

if(ADIOS2_RUN_MPI_MPMD_TESTS)
  if(ADIOS2_HAVE_SST)
    include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSTBPMxM.cmake)
    include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSTBPMxN.cmake)
    include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSTBPMx1.cmake)
    include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSTFFSMxM.cmake)
    include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSTFFSMxN.cmake)
    include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSTFFSMx1.cmake)
    if(ADIOS2_SST_HAVE_LIBFABRIC)
    #  include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSTBPRDMAMxN.cmake)
    #  include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSTFFSRDMAMxN.cmake)
    endif()
  endif()
  if(NOT MSVC)
    include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSCMxM.cmake)
    include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSCMxN.cmake)
    include(${CMAKE_CURRENT_SOURCE_DIR}/TestSSCMx1.cmake)
  endif()
endif()
