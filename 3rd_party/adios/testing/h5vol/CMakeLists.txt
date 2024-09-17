#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if(ADIOS2_HAVE_MPI AND HDF5_IS_PARALLEL)
  set(hdf5_mpi MPI_ONLY)
  add_definitions(-DTEST_HDF5_MPI)
  set(hdf5_sfx ".MPI")
else()
  set(hdf5_mpi MPI_NONE)
  set(hdf5_sfx ".Serial")
endif()



gtest_add_tests_helper(WriteReadBPFile ${hdf5_mpi} H5Vol Util.H5VOL. "")
#target_link_libraries(Test.Util.H5VOL.WriteReadBPFile$  ${HDF5_C_LIBRARIES})


if(HDF5_C_INCLUDE_DIRS)
  target_include_directories(Test.Util.H5VOL.WriteReadBPFile${hdf5_sfx}
    PRIVATE ${HDF5_C_INCLUDE_DIRS}
  )
else()
  target_include_directories(Test.Util.H5VOL.WriteReadBPFile${hdf5_sfx}
    PRIVATE ${HDF5_INCLUDE_DIRS}
  )
endif()
target_link_libraries(Test.Util.H5VOL.WriteReadBPFile${hdf5_sfx} ${HDF5_C_LIBRARIES} adios2_h5vol)

