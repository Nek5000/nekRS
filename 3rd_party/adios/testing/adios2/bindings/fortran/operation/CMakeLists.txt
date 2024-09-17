if(ADIOS2_HAVE_SZ)
  fortran_add_test_helper(BPWriteReadSZ2D MPI_ONLY)
  fortran_add_test_helper(BPWriteReadSZ3D MPI_ONLY)
endif()

if(ADIOS2_HAVE_ZFP)
  fortran_add_test_helper(BPWriteReadZfp2D MPI_ONLY)
  fortran_add_test_helper(BPWriteReadZfp2DRemove MPI_ONLY)
endif()
