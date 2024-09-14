#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

include(ADIOSFunctions)

if(ADIOS2_HAVE_MPI AND NOT MSVC)
  gtest_add_tests_helper(Base MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscBase.MPI "" TRUE)

  gtest_add_tests_helper(3DMemSelect MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSsc3DMemSelect.MPI "" TRUE)

  gtest_add_tests_helper(Struct MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscStruct.MPI "" TRUE)

  gtest_add_tests_helper(SingleStep MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscSingleStep.MPI "" TRUE)

  gtest_add_tests_helper(SuperLarge MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscSuperLarge.MPI "" TRUE)

  gtest_add_tests_helper(String MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscString.MPI "" TRUE)

  gtest_add_tests_helper(BaseUnlocked MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscBaseUnlocked.MPI "" TRUE)

  gtest_add_tests_helper(LockBeforeEndStep MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscLockBeforeEndStep.MPI "" TRUE)

  gtest_add_tests_helper(OnlyOneStep MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscOnlyOneStep.MPI "" TRUE)

  gtest_add_tests_helper(OnlyTwoSteps MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscOnlyTwoSteps.MPI "" TRUE)

  gtest_add_tests_helper(Unbalanced MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscUnbalanced.MPI "" TRUE)

  gtest_add_tests_helper(WriterMultiblock MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscWriterMultiblock.MPI "" TRUE)

  gtest_add_tests_helper(ReaderMultiblock MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscReaderMultiblock.MPI "" TRUE)

  gtest_add_tests_helper(NoAttributes MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscNoAttributes.MPI "" TRUE)

  gtest_add_tests_helper(NoSelection MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscNoSelection.MPI "" TRUE)

  gtest_add_tests_helper(7d MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSsc7d.MPI "" TRUE)

  gtest_add_tests_helper(MoreReadersThanWriters MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscMoreReadersThanWriters.MPI "" TRUE)

  gtest_add_tests_helper(MoreWritersThanReaders MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscMoreWritersThanReaders.MPI "" TRUE)

  gtest_add_tests_helper(Xgc2Way MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscXgc2Way.MPI "" TRUE)

  gtest_add_tests_helper(Xgc3Way MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscXgc3Way.MPI "" TRUE)

  gtest_add_tests_helper(Xgc3WayMatchedSteps MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscXgc3WayMatchedSteps.MPI "" TRUE)

  gtest_add_tests_helper(VaryingSteps MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscVaryingSteps.MPI "" TRUE)

  gtest_add_tests_helper(ZeroBlock MPI_ONLY Ssc Engine.SSC. "")
  SetupTestPipeline(Engine.SSC.SscEngineTest.TestSscZeroBlock.MPI "" TRUE)

endif()
