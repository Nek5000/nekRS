program TestBPFortranToCppWriter
   use mpi
   use adios2
   implicit none

#ifndef __GFORTRAN__
#if !defined(__GNUC__) || defined(_CRAYFTN)
   interface
      integer function iargc()
      end function iargc
   end interface
#endif
#endif

   ! command line input
   character(len=256) :: engine
   integer :: numargs

   integer, parameter :: inx = 10

   integer(kind=8), dimension(1) :: count_dims, changing_count_dims
   integer(kind=8), dimension(2) :: count2_dims
   integer :: irank, isize, ierr, i, s

   type(adios2_adios) :: adios
   type(adios2_io) :: ioWrite
   type(adios2_variable) :: vGlobalValue, vLA_1D, vLA_2D, vLA_1D_changing
   type(adios2_engine) :: bpWriter
   integer(kind=8) current_step

   ! local arrays
   integer(kind=8), dimension(inx) :: LA_1D
   integer(kind=8), dimension(inx) :: LA_1D_changing
   integer(kind=8), dimension(inx,2) :: LA_2D

   character(len=:), allocatable :: engine_type

   ! Get command line input
   numargs = iargc()
   ! print *,"Number of arguments:",numargs
   if ( numargs > 0 ) then
      call getarg(1, engine)
   endif

   ! Launch MPI
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

   ! Variable dimensions
   count_dims(1) = inx
   changing_count_dims(1) = inx
   count2_dims(1) = inx
   count2_dims(2) = 2

   call adios2_init(adios, MPI_COMM_WORLD, ierr)
   call adios2_declare_io(ioWrite, adios, "ioWrite", ierr)

   call adios2_define_variable(vGlobalValue, ioWrite, "inx", &
      adios2_type_integer4,  ierr)
   call adios2_define_variable(vLA_1D, ioWrite, "localarray_1D", &
      adios2_type_integer8, 1, &
      adios2_null_dims, adios2_null_dims, count_dims, &
      .false., ierr)
   call adios2_define_variable(vLA_2D, ioWrite, "localarray_2D", &
      adios2_type_integer8, 2, &
      adios2_null_dims, adios2_null_dims, count2_dims, &
      adios2_variable_dims, ierr)
   call adios2_define_variable(vLA_1D_changing, ioWrite, "localarray_1D_changing", &
      adios2_type_integer8, 1, &
      adios2_null_dims, adios2_null_dims, changing_count_dims, &
      adios2_variable_dims, ierr)

   call adios2_set_engine(ioWrite, trim(engine), ierr)
   call adios2_io_engine_type(engine_type, ioWrite, ierr)
   if (irank == 0) print *,"engine type :",trim(engine_type)


#if ADIOS2_USE_MPI
   call adios2_open(bpWriter, ioWrite, "FortranToCpp_MPI.bp", &
      adios2_mode_write, ierr)
#else
   call adios2_open(bpWriter, ioWrite, "FortranToCpp.bp", &
      adios2_mode_write, ierr)
#endif

   do s = 1, 3
      call adios2_begin_step(bpWriter, ierr)

      call adios2_current_step(current_step, bpWriter, ierr)
      if (current_step /= s - 1) then
         write(*,*) 'wrong current step'
         stop 1
      end if

      if (irank == 0 .and. s == 1) then
         call adios2_put(bpWriter, vGlobalValue, inx, ierr)
      end if

      do i = 1, inx
         LA_1D(i) = (irank+1)*1000 + s*100 + i
         LA_2D(i,1) = (irank+1)*1000 + s*100 + i
         LA_2D(i,2) = LA_2D(i,1) + inx
      end do

      ! changing count
      changing_count_dims(1) = s
      call adios2_set_selection(vLA_1D_changing, 1, adios2_null_dims, &
         changing_count_dims, ierr)
      do i = 1, s
         LA_1D_changing(i) = LA_1D(i) + s
      end do

      call adios2_put(bpWriter, vLA_1D, LA_1D, ierr)
      call adios2_put(bpWriter, vLA_2D, LA_2D, ierr)
      call adios2_put(bpWriter, vLA_1D_changing, LA_1D_changing, ierr)

      call adios2_end_step(bpWriter, ierr)
   end do

   call adios2_close(bpWriter, ierr)
   call adios2_finalize(adios, ierr)
   call MPI_Finalize(ierr)

end program TestBPFortranToCppWriter

