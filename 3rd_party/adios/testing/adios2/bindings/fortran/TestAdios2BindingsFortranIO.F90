
! ======================================================================
module testing_adios
! ======================================================================
  use adios2
  use mpi
  implicit none

  type(adios2_adios) :: adios
contains
  ! ------------------------------------------------------------
  subroutine testing_adios_init()
  ! ------------------------------------------------------------

    integer :: ierr

    call adios2_init(adios, MPI_COMM_WORLD, ierr)

  end subroutine testing_adios_init

  ! ------------------------------------------------------------
  subroutine testing_adios_finalize()
  ! ------------------------------------------------------------

    integer :: ierr

    call adios2_finalize(adios, ierr)

  end subroutine testing_adios_finalize

end module testing_adios

! ======================================================================
module testing_adios_io
! ======================================================================
  use testing_adios
  use adios2
  implicit none

  type(adios2_io) :: io

contains
  ! ------------------------------------------------------------
  subroutine testing_adios_io_init()
  ! ------------------------------------------------------------

    integer :: ierr

    call testing_adios_init()
    call adios2_declare_io(io, adios, "TestIo", ierr)

  end subroutine testing_adios_io_init

  ! ------------------------------------------------------------
  subroutine testing_adios_io_finalize()
  ! ------------------------------------------------------------

    integer :: ierr
    logical :: result

    ! FIXME, shouldn't we be able to do this by handle?
    call adios2_remove_io(result, adios, "TestIo", ierr)
    if ((ierr /= 0) .or. (result .neqv. .true.)) then
       write(*,*) "FAIL: adios2_remove_io"
       stop 1
    end if
    call testing_adios_finalize()

  end subroutine testing_adios_io_finalize

end module testing_adios_io

! ======================================================================

! ------------------------------------------------------------
subroutine testing_adios_io_engine()
! ------------------------------------------------------------
! test engine related functionality that's part of IO

  use testing_adios_io
  implicit none

  integer :: ierr
  type(adios2_engine) :: engine

  character(len=:), allocatable :: engine_type

  call testing_adios_io_init()

  ! Engine related functionality
  call adios2_set_engine(io, "file", ierr)

  call adios2_io_engine_type(engine_type, io, ierr)
  if (engine_type /= "file") then
     write(*,*) "FAIL adios2_io_engine_type"
     stop 1
  end if
  deallocate(engine_type)

  call adios2_open(engine, io, "ftypes.bp", adios2_mode_write, ierr)

  if (engine%type /= "BP5Writer") then
     write(*,*) "FAIL engine%type"
     stop 1
  end if
  ! // FIXME, I'd like to check that the engine type itself is correct, but
  ! // there's no (function-style) API to get it
  ! // FIXME, I'd like to check the engine's name, but there's no API to get it

  call adios2_io_engine_type(engine_type, io, ierr)
  if (engine_type /= "file") then
     write(*,*) "FAIL adios2_io_engine_type"
     stop 1
  end if
  deallocate(engine_type)

  call testing_adios_io_finalize()

end subroutine testing_adios_io_engine

! ------------------------------------------------------------
subroutine testing_adios_io_engine_default()
! ------------------------------------------------------------
  use testing_adios_io
  implicit none

  integer :: ierr
  type(adios2_engine) :: engine

  character(len=:), allocatable :: engine_type

  call testing_adios_io_init()

  ! Engine related functionality
  call adios2_set_engine(io, "", ierr)

  call adios2_io_engine_type(engine_type, io, ierr)
  if (engine_type /= "") then
     write(*,*) "FAIL adios2_io_engine_type"
     stop 1
  end if
  deallocate(engine_type)

  call adios2_open(engine, io, "ftypes.bp", adios2_mode_write, ierr)

  if (engine%type /= "BP5Writer") then
     write(*,*) "FAIL engine%type"
     stop 1
  end if
  ! // FIXME, I'd like to check that the engine type itself is correct, but
  ! // there's no (function-style) API to get it
  ! // FIXME, I'd like to check the engine's name, but there's no API to get it

  call adios2_io_engine_type(engine_type, io, ierr)
  if (engine_type /= "") then
     write(*,*) "FAIL adios2_io_engine_type"
     stop 1
  end if
  deallocate(engine_type)

  call testing_adios_io_finalize

end subroutine testing_adios_io_engine_default

! ======================================================================
program main
! ======================================================================
  use mpi
  implicit none

  integer :: ierr
  external testing_adios_io_engine
  external testing_adios_io_engine_default

  INTEGER provided

  ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
  call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)

  call testing_adios_io_engine()
  call testing_adios_io_engine_default()

  call MPI_Finalize(ierr)

end program main
