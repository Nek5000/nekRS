program adios_fortran_test
  use adios2
  integer :: ierr, irank, isize
  type(adios2_adios) :: adios

  call adios2_init(adios, ierr)
  call adios2_finalize(adios, ierr)
end program adios_fortran_test
