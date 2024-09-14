program adios_fortran_test
  use mpi
  use adios2
  integer :: ierr, irank, isize
  type(adios2_adios) :: adios

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)

  call adios2_init(adios, MPI_COMM_WORLD, ierr)
  call adios2_finalize(adios, ierr)

  call MPI_Finalize(ierr)
end program adios_fortran_test
