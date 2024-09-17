!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_adios_init_mod.F90 : ADIOS2 Fortran bindings for ADIOS class Init
!                              functions
!

module adios2_adios_init_mod
#ifdef ADIOS2_HAVE_FORTRAN_SUBMODULES
    use adios2_parameters_mod
    implicit none

    interface adios2_init
        module subroutine adios2_init_serial(adios, ierr)
            type(adios2_adios), intent(out) :: adios
            integer, intent(out) :: ierr
        end subroutine
        module subroutine adios2_init_config_serial(adios, config_file, ierr)
            type(adios2_adios), intent(out) :: adios
            character*(*), intent(in) :: config_file
            integer, intent(out) :: ierr
        end subroutine

#ifdef ADIOS2_HAVE_MPI_F
        module subroutine adios2_init_mpi(adios, comm, ierr)
            type(adios2_adios), intent(out) :: adios
            integer, intent(in) :: comm
            integer, intent(out) :: ierr
        end subroutine
        module subroutine adios2_init_config_mpi(adios, config_file, comm, ierr)
            type(adios2_adios), intent(out) :: adios
            character*(*), intent(in) :: config_file
            integer, intent(in) :: comm
            integer, intent(out) :: ierr
        end subroutine
#endif
    end interface
#else

    use adios2_adios_init_serial_mod

#ifdef ADIOS2_HAVE_MPI_F
    use adios2_adios_init_mpi_mod
#endif /*ADIOS2_HAVE_MPI_F*/
#endif /*ADIOS2_HAVE_FORTRAN_SUBMODULES*/

end module
