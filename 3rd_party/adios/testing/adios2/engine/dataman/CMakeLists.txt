#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

foreach(tst IN ITEMS
  1D
  1DSuperLarge
  2DMemSelect 3DMemSelect
  WriterDoubleBuffer WriterSingleBuffer
  ReaderDoubleBuffer ReaderSingleBuffer
  Reliable
  )
  gtest_add_tests_helper(${tst} MPI_NONE DataMan Engine.DataMan. "")
  set_tests_properties(${Test.Engine.DataMan.${tst}-TESTS}
    PROPERTIES RUN_SERIAL TRUE
  )
endforeach()

if (ADIOS2_HAVE_Python)
    python_add_test(NAME Test.Engine.DataMan1D.Serial SCRIPT TestDataMan1D.py)
    python_add_test(NAME Test.Engine.DataMan1xN.Serial SCRIPT TestDataMan1xN.py)
    python_add_test(NAME Test.Engine.DataManSingleValues SCRIPT TestDataManSingleValues.py)
endif()

if(ADIOS2_HAVE_ZFP)
    gtest_add_tests_helper(2DZfp MPI_NONE DataMan Engine.DataMan. "")
    set_tests_properties(${Test.Engine.DataMan.2DZfp-TESTS}
        PROPERTIES RUN_SERIAL TRUE
        )
endif()

if(ADIOS2_HAVE_SZ)
    gtest_add_tests_helper(2DSz MPI_NONE DataMan Engine.DataMan. "")
    set_tests_properties(${Test.Engine.DataMan.2DSz-TESTS}
        PROPERTIES RUN_SERIAL TRUE
        )
endif()
