*************
Running Tests
*************

ADIOS2 uses `googletest <https://github.com/google/googletest>`_ to enable automatic testing after a CMake build. To run tests just type after building with make, run:

.. code-block:: bash

    $ ctest    
      or  
    $ make test

The following screen will appear providing information on the status of each finalized test:

.. code-block:: bash

    Test project /home/wfg/workspace/build
    Start  1: ADIOSInterfaceWriteTest.DefineVarChar1x10
    1/46 Test  #1: ADIOSInterfaceWriteTest.DefineVarChar1x10 .......................   Passed    0.06 sec
    Start  2: ADIOSInterfaceWriteTest.DefineVarShort1x10
    2/46 Test  #2: ADIOSInterfaceWriteTest.DefineVarShort1x10 ......................   Passed    0.04 sec
    Start  3: ADIOSInterfaceWriteTest.DefineVarInt1x10
    3/46 Test  #3: ADIOSInterfaceWriteTest.DefineVarInt1x10 ........................   Passed    0.04 sec
    Start  4: ADIOSInterfaceWriteTest.DefineVarLong1x10
    ... 
    128/130 Test #128: ADIOSZfpWrapper.UnsupportedCall ..........................................   Passed    0.05 sec
        Start 129: ADIOSZfpWrapper.MissingMandatoryParameter
    129/130 Test #129: ADIOSZfpWrapper.MissingMandatoryParameter ................................   Passed    0.05 sec
        Start 130: */TestManyVars.DontRedefineVars/*
    130/130 Test #130: */TestManyVars.DontRedefineVars/* ........................................   Passed    0.08 sec

    100% tests passed, 0 tests failed out of 130

    Total Test time (real) = 204.82 sec



