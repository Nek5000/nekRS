#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

# These tests are longer than others but should still be pretty fast
set(CTEST_TEST_TIMEOUT 60)

set(BP3_DIR ${CMAKE_CURRENT_BINARY_DIR}/bp3)
set(BP4_DIR ${CMAKE_CURRENT_BINARY_DIR}/bp4)
set(BP5_DIR ${CMAKE_CURRENT_BINARY_DIR}/bp5)
file(MAKE_DIRECTORY ${BP3_DIR})
file(MAKE_DIRECTORY ${BP4_DIR})
file(MAKE_DIRECTORY ${BP5_DIR})

if(ADIOS2_HAVE_SZ)
  bp_gtest_add_tests_helper(WriteReadSZ MPI_ALLOW)
  gtest_add_tests_helper(SzComplex MPI_ALLOW BPWriteRead Engine. "")
endif()


if(ADIOS2_HAVE_ZFP)
  bp_gtest_add_tests_helper(WriteReadZfp MPI_ALLOW)
  bp_gtest_add_tests_helper(WriteReadZfpHighLevelAPI MPI_ALLOW)
  bp3_bp4_gtest_add_tests_helper(WriteReadZfpRemoveOperations MPI_ALLOW)
  bp_gtest_add_tests_helper(WriteReadZfpConfig MPI_ALLOW)

  gtest_add_tests_helper(ZfpComplex MPI_ALLOW BPWriteRead Engine. "")

  foreach(tgt
      ${Test.Engine.BP.WriteReadZfpConfig-TARGETS}
      )
    target_compile_definitions(${tgt} PRIVATE
      "XML_CONFIG_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
      )
  endforeach()

  if(ADIOS2_HAVE_CUDA AND ADIOS2_HAVE_ZFP_CUDA)
    gtest_add_tests_helper(WriteReadZfpCuda MPI_ALLOW BP Engine.BP. .BP4
      WORKING_DIRECTORY ${BP4_DIR} EXTRA_ARGS "BP4"
      )

    foreach(tgt ${Test.Engine.BP.WriteReadZfpCuda-TARGETS})
      target_sources(${tgt} PRIVATE CudaRoutines.cu)
      target_link_libraries(${tgt} CUDA::cudart)
    endforeach()
  endif()
endif()

if(ADIOS2_HAVE_MGARD)
  bp_gtest_add_tests_helper(WriteReadMGARD MPI_ALLOW)
  bp_gtest_add_tests_helper(WriteReadMGARDPlus MPI_ALLOW)
  if(ADIOS2_HAVE_CUDA AND MGARD_ENABLE_CUDA)
    gtest_add_tests_helper(WriteReadMGARDCuda MPI_ALLOW BP Engine.BP. .BP4
      WORKING_DIRECTORY ${BP4_DIR} EXTRA_ARGS "BP4"
      )

    foreach(tgt ${Test.Engine.BP.WriteReadMGARDCuda-TARGETS})
      target_sources(${tgt} PRIVATE CudaRoutines.cu)
      target_link_libraries(${tgt} CUDA::cudart)
    endforeach()
  endif()

  if(ADIOS2_HAVE_MGARD_MDR)
    gtest_add_tests_helper(WriteReadMGARDMDR MPI_ALLOW BP Engine.BP. .BP5
    WORKING_DIRECTORY ${BP5_DIR} EXTRA_ARGS "BP5"
    )
  endif()
endif()

if(ADIOS2_HAVE_BZip2)
  bp_gtest_add_tests_helper(WriteReadBZIP2 MPI_ALLOW)
endif()

if(ADIOS2_HAVE_PNG)
  bp_gtest_add_tests_helper(WriteReadPNG MPI_ALLOW)
endif()

if(ADIOS2_HAVE_Blosc2)
  bp_gtest_add_tests_helper(WriteReadBlosc2 MPI_ALLOW)
endif()
