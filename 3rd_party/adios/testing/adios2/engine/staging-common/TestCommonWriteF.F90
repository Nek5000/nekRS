subroutine usage()
    print *, "Usage: TestCommonWriteF engine filename"
end subroutine usage


program TestSstWrite
  use sst_test_data
#if ADIOS2_USE_MPI
  use mpi
#endif
  use adios2
  implicit none
  external usage

#if defined(ADIOS2_HAVE_FORTRAN_F03_ARGS)
# define ADIOS2_ARGC() command_argument_count()
# define ADIOS2_ARGV(i, v) call get_command_argument(i, v)
#elif defined(ADIOS2_HAVE_FORTRAN_GNU_ARGS)
# define ADIOS2_ARGC() iargc()
# define ADIOS2_ARGV(i, v) call getarg(i, v)
#else
# define ADIOS2_ARGC() 1
# define ADIOS2_ARGV(i, v)
#endif

  integer :: numargs

  integer(kind = 8), dimension(1)::shape_dims, start_dims, count_dims
  integer(kind = 8), dimension(2)::shape_dims2, start_dims2, count_dims2
  integer(kind = 8), dimension(2)::shape_dims3, start_dims3, count_dims3
  integer(kind = 8), dimension(1)::shape_time, start_time, count_time
  integer::irank, isize, ierr, i, insteps, status

  character(len=256) :: filename, engine, params

  type(adios2_adios)::adios
  type(adios2_io)::ioWrite
  type(adios2_variable), dimension(:), allocatable :: variables
  type(adios2_engine)::sstWriter;

  !read handlers
  integer(kind = 8)::localtime

#if ADIOS2_USE_MPI
  integer::testComm, color, key, provided, threadSupportLevel
#endif

  allocate(variables(20))

  numargs = ADIOS2_ARGC()
  if ( numargs < 2 ) then
     call usage()
     call exit(1)
  endif


  ADIOS2_ARGV(1, engine)
  ADIOS2_ARGV(2, filename)
  if ( numargs > 2 ) then
     ADIOS2_ARGV(3, params)
  endif

#if ADIOS2_USE_MPI
  threadSupportLevel = MPI_THREAD_SINGLE;
  if (engine == "SST") then
      threadSupportLevel = MPI_THREAD_MULTIPLE;
  endif

  !Launch MPI

  ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
  call MPI_Init_thread(threadSupportLevel, provided, ierr)

  call MPI_Comm_rank(MPI_COMM_WORLD, key, ierr);

  color = 1
  call MPI_Comm_split(MPI_COMM_WORLD, color, key, testComm, ierr);

  call MPI_Comm_rank(testComm, irank, ierr)
  call MPI_Comm_size(testComm, isize, ierr)
#else
  ! No MPI
  irank = 0;
  isize = 1;
#endif

  !Application variables
  insteps = 10;

  !Variable dimensions
  shape_dims(1) = isize * nx
  start_dims(1) = irank * nx
  count_dims(1) = nx

  shape_dims2 = (/ 2, isize *nx /)
  start_dims2 = (/ 0, irank *nx /)
  count_dims2 = (/ 2, nx /)

  shape_dims3 = (/ isize *nx, 2 /)
  start_dims3 = (/ irank *nx, 0 /)
  count_dims3 = (/ nx, 2 /)

  shape_time = (/ isize /)
  start_time = (/ irank /)
  count_time = (/ 1 /)

#if ADIOS2_USE_MPI
  !Create adios handler passing the communicator and error flag
  call adios2_init(adios, testComm, ierr)
#else
  call adios2_init(adios, ierr)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!WRITER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!Declare an IO process configuration inside adios
  call adios2_declare_io(ioWrite, adios, "ioWrite", ierr)

  if (numargs > 2) then
     call adios2_set_parameters(ioWrite, params, ierr)
  endif

  call adios2_set_engine(ioWrite, engine, ierr)

  !Defines a variable to be written
  call adios2_define_variable(variables(12), ioWrite, "scalar_r64", &
       adios2_type_dp, ierr)

  call adios2_define_variable(variables(1), ioWrite, "i8", &
       adios2_type_integer1, 1, &
       shape_dims, start_dims, count_dims, &
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(2), ioWrite, "i16", &
       adios2_type_integer2, 1, &
       shape_dims, start_dims, count_dims, &
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(3), ioWrite, "i32", &
       adios2_type_integer4, 1, &
       shape_dims, start_dims, count_dims, &
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(4), ioWrite, "i64", &
       adios2_type_integer8, 1, &
       shape_dims, start_dims, count_dims, &
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(5), ioWrite, "r32", &
       adios2_type_real, 1, &
       shape_dims, start_dims, count_dims,&
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(6), ioWrite, "r64", &
       adios2_type_dp, 1, &
       shape_dims, start_dims, count_dims, &
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(7), ioWrite, "r64_2d", &
       adios2_type_dp, 2, &
       shape_dims2, start_dims2, count_dims2, &
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(8), ioWrite, "r64_2d_rev", &
       adios2_type_dp, 2, &
       shape_dims3, start_dims3, count_dims3, &
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(9), ioWrite, "time", &
       adios2_type_integer8, 1, &
       shape_time, start_time, count_time, &
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(10), ioWrite, "c32", &
       adios2_type_complex, 1, &
       shape_dims, start_dims, count_dims,&
       adios2_constant_dims, ierr)

  call adios2_define_variable(variables(11), ioWrite, "c64", &
       adios2_type_complex_dp, 1, &
       shape_dims, start_dims, count_dims, &
       adios2_constant_dims, ierr)

  call adios2_open(sstWriter, ioWrite, filename, adios2_mode_write, ierr)

  !Put array contents to bp buffer, based on var1 metadata
  do i = 1, insteps
     call GenerateTestData(i - 1, irank)
     call adios2_begin_step(sstWriter, adios2_step_mode_append, -1.0, &
                            status, ierr)
     call adios2_put(sstWriter, variables(12), data_scalar_r64, ierr)
     call adios2_put(sstWriter, variables(1), data_I8, ierr)
     call adios2_put(sstWriter, variables(2), data_I16, ierr)
     call adios2_put(sstWriter, variables(3), data_I32, ierr)
     call adios2_put(sstWriter, variables(4), data_I64, ierr)
     call adios2_put(sstWriter, variables(5), data_R32, ierr)
     call adios2_put(sstWriter, variables(6), data_R64, ierr)
     call adios2_put(sstWriter, variables(7), data_R64_2d, ierr)
     call adios2_put(sstWriter, variables(8), data_R64_2d_rev, ierr)
#if !defined(_CRAYFTN)
     localtime = 0    ! should be time(), but non-portable and value is unused
     call adios2_put(sstWriter, variables(9), loc(localtime), ierr)
#endif
     call adios2_put(sstWriter, variables(10), data_C32, ierr)
     call adios2_put(sstWriter, variables(11), data_C64, ierr)
     call adios2_end_step(sstWriter, ierr)
  end do

  !Closes engine1 and deallocates it, becomes unreachable
  call adios2_close(sstWriter, ierr)

  !Deallocates adios and calls its destructor
  call adios2_finalize(adios, ierr)

#if ADIOS2_USE_MPI
#ifdef CRAY_MPICH_VERSION
  call MPI_Barrier(MPI_COMM_WORLD)
#else
  call MPI_Finalize(ierr)
#endif
#endif

 end program TestSstWrite
