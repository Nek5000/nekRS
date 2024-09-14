#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

include(ADIOSFunctions)

add_executable(Test.Engine.Skeleton.Writer TestSkeletonWriter.cpp)

add_executable(Test.Engine.Skeleton.Reader TestSkeletonReader.cpp)

if(ADIOS2_HAVE_MPI)
  target_link_libraries(Test.Engine.Skeleton.Writer adios2::cxx11_mpi MPI::MPI_C)
  target_link_libraries(Test.Engine.Skeleton.Reader adios2::cxx11_mpi MPI::MPI_C)
  set(cmd_executor
    ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS} ${MPIEXEC_NUMPROC_FLAG} 1
  )
else()
  target_link_libraries(Test.Engine.Skeleton.Writer adios2::cxx11)
  target_link_libraries(Test.Engine.Skeleton.Reader adios2::cxx11)
  set(cmd_executor)
endif()

add_test(NAME Engine.Skeleton.Writer
  COMMAND ${cmd_executor} $<TARGET_FILE:Test.Engine.Skeleton.Writer> 
)

add_test(NAME Engine.Skeleton.Reader
  COMMAND ${cmd_executor} $<TARGET_FILE:Test.Engine.Skeleton.Reader>
)
SetupTestPipeline(Engine.Skeleton "Writer;Reader" TRUE)

add_test(NAME Engine.Skeleton.Writer.Validate
  COMMAND ${DIFF_COMMAND} -u -w
    ${CMAKE_CURRENT_SOURCE_DIR}/TestSkeletonWriterExpectedOutput.txt
    TestSkeletonWriterOutput.txt
)
SetupTestPipeline(Engine.Skeleton "Writer;Writer.Validate" FALSE)

add_test(NAME Engine.Skeleton.Reader.Validate
  COMMAND ${DIFF_COMMAND} -u -w
    ${CMAKE_CURRENT_SOURCE_DIR}/TestSkeletonReaderExpectedOutput.txt
    TestSkeletonReaderOutput.txt
)
SetupTestPipeline(Engine.Skeleton "Reader;Reader.Validate" FALSE)
