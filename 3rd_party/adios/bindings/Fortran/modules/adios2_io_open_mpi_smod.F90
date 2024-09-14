!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_io_open_mod_mpi.F90 : ADIOS2 Fortran bindings for IO
!                               class open function (MPI variants)
!

#ifdef ADIOS2_HAVE_FORTRAN_SUBMODULES
# define ADIOS2_MODULE_PROCEDURE module
#else
# define ADIOS2_MODULE_PROCEDURE
#endif

#ifdef ADIOS2_HAVE_FORTRAN_SUBMODULES
submodule ( adios2_io_open_mod ) adios2_io_open_mpi_smod
#else
module adios2_io_open_mpi_mod
#endif

    use adios2_parameters_mod
    implicit none

#ifndef ADIOS2_HAVE_FORTRAN_SUBMODULES
    interface adios2_open
        module procedure adios2_open_new_comm
    end interface
#endif
    external adios2_open_new_comm_f2c
    external adios2_engine_get_type_f2c

contains

    ADIOS2_MODULE_PROCEDURE subroutine adios2_open_new_comm( &
            engine, io, name, adios2_mode, comm, ierr)
        type(adios2_engine), intent(out) :: engine
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer, intent(in) :: adios2_mode
        integer, intent(in) :: comm
        integer, intent(out) :: ierr

        call adios2_open_new_comm_f2c(engine%f2c, io%f2c, &
                                      TRIM(ADJUSTL(name))//char(0), &
                                      adios2_mode, comm, ierr)

        if( ierr == 0 ) then
            engine%valid = .true.
            engine%name = name
            call adios2_engine_get_type_f2c(engine%type, engine%f2c, ierr)
            engine%mode = adios2_mode
        end if

    end subroutine

#ifdef ADIOS2_HAVE_FORTRAN_SUBMODULES
end submodule
#else
end module
#endif
