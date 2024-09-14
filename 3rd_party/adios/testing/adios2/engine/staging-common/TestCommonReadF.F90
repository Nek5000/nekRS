!!***************************
subroutine usage()
    print *, "Usage: TestCommonReadF engine filename"
end subroutine usage


program TestSstRead
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

  integer(kind = 8), dimension(1)::start_dims, count_dims
  integer(kind = 8), dimension(2)::start_dims2, count_dims2
  integer(kind = 8), dimension(2)::start_dims3, count_dims3
  integer:: irank, isize, ierr, insteps, step_status

  character(len=256) :: filename, engine, params

  integer :: writerSize, myStart, myLength

  type(adios2_adios)::adios
  type(adios2_io)::ioRead
  type(adios2_variable), dimension(:), allocatable :: variables
  type(adios2_engine)::sstReader;

  !read handlers
  integer::ndims
  integer(kind = 8), dimension(:), allocatable::shape_in

#if ADIOS2_USE_MPI
  integer::key, color, testComm, provided, threadSupportLevel
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

  insteps = 0;

#if ADIOS2_USE_MPI
  threadSupportLevel = MPI_THREAD_SINGLE;
  if (engine == "SST") then
      threadSupportLevel = MPI_THREAD_MULTIPLE;
  endif

  !Launch MPI

  ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
  call MPI_Init_thread(threadSupportLevel, provided, ierr)

  call MPI_Comm_rank(MPI_COMM_WORLD, key, ierr);

  color = 2
  call MPI_Comm_split(MPI_COMM_WORLD, color, key, testComm, ierr);

  call MPI_Comm_rank(testComm, irank, ierr)
  call MPI_Comm_size(testComm, isize, ierr)

  !Create adios handler passing the communicator and error flag
  call adios2_init(adios, testComm, ierr)
#else
  irank = 0;
  isize = 1;

  !Create adios handler passing the error flag
  call adios2_init(adios, ierr)
#endif
  !!!!!!!!!!!!!!!!!!!!!!!!! READER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Declare io reader
  call adios2_declare_io(ioRead, adios, "ioRead", ierr)

  if (numargs > 2) then
     call adios2_set_parameters(ioRead, params, ierr)
  endif
  call adios2_set_engine(ioRead, engine, ierr)

  ! Open sstReader engine
  call adios2_open(sstReader, ioRead, filename, adios2_mode_read, ierr)

  do
    call adios2_begin_step(sstReader, adios2_step_mode_read, -1., step_status, ierr)

    if(step_status == adios2_step_status_end_of_stream) exit


    call adios2_inquire_variable(variables(1), ioRead, "i8", ierr)
    if (variables(1)%name /= 'i8') then
       write(*,*) 'i8 not recognized'
       stop 1
    end if
    if (variables(1)%type /= adios2_type_integer1) then
       write(*,*) 'i8 type not recognized'
       stop 1
    end if
    call adios2_variable_shape(shape_in, ndims, variables(1),  ierr)
    if (ndims /= 1) then
       write(*,*) 'i8 ndims is not 1'
       stop 1
    end if
    if (modulo(shape_in(1), int(nx, 8)) /= 0) then
       write(*,*) 'i8 shape_in read failed'
       stop 1
    end if
    writerSize = INT(shape_in(1)) / nx
    deallocate(shape_in)

    call adios2_inquire_variable(variables(2), ioRead, "i16", ierr)
    if (variables(2)%name /= 'i16') then
       write(*,*) 'i16 not recognized'
       stop 1
    end if
    if (variables(2)%type /= adios2_type_integer2) then
       write(*,*) 'i16 type not recognized'
       stop 1
    end if
    call adios2_variable_shape( shape_in, ndims, variables(2), ierr)
    if (ndims /= 1) then
       write(*,*) 'i16 ndims is not 1'
       stop 1
    end if
    if (shape_in(1) /= nx*writerSize) then
       write(*,*) 'i16 shape_in read failed'
       stop 1
    end if
    deallocate(shape_in)

    call adios2_inquire_variable(variables(3), ioRead, "i32", ierr)
    if (variables(3)%name /= 'i32') then
       write(*,*) 'i32 not recognized'
       stop 1
    end if
    if (variables(3)%type /= adios2_type_integer4) then
       write(*,*) 'i32 type not recognized'
       stop 1
    end if
    call adios2_variable_shape(shape_in, ndims, variables(3), ierr)
    if (ndims /= 1) then
       write(*,*) 'i32 ndims is not 1'
       stop 1
    end if
    if (shape_in(1) /= nx*writerSize) then
       write(*,*) 'i32 shape_in read failed'
       stop 1
    end if
    deallocate(shape_in)

    call adios2_inquire_variable(variables(4), ioRead, "i64", ierr)
    if (variables(4)%name /= 'i64') then
       write(*,*) 'i64 not recognized'
       stop 1
    end if
    if (variables(4)%type /= adios2_type_integer8) then
       write(*,*) 'i64 type not recognized'
       stop 1
    end if
    call adios2_variable_shape(shape_in, ndims, variables(4), ierr)
    if (ndims /= 1) then
       write(*,*) 'i64 ndims is not 1'
       stop 1
    end if
    if (shape_in(1) /= nx*writerSize) then
       write(*,*) 'i64 shape_in read failed'
       stop 1
    end if
    deallocate(shape_in)

    call adios2_inquire_variable(variables(5), ioRead, "r32", ierr)
    if (variables(5)%name /= 'r32') then
       write(*,*) 'r32 not recognized'
       stop 1
    end if
    if (variables(5)%type /= adios2_type_real) then
       write(*,*) 'r32 type not recognized'
       stop 1
    end if
    call adios2_variable_shape(shape_in, ndims, variables(5), ierr)
    if (ndims /= 1) then
       write(*,*) 'r32 ndims is not 1'
       stop 1
    end if
    if (shape_in(1) /= nx*writerSize) then
       write(*,*) 'r32 shape_in read failed'
       stop 1
    end if
    deallocate(shape_in)

    call adios2_inquire_variable(variables(6), ioRead, "r64", ierr)
    if (variables(6)%name /= 'r64') then
       write(*,*) 'r64 not recognized'
       stop 1
    end if
    if (variables(6)%type /= adios2_type_dp) then
       write(*,*) 'r64 type not recognized'
       stop 1
    end if
    call adios2_variable_shape(shape_in, ndims, variables(6), ierr)
    if (ndims /= 1) then
       write(*,*) 'r64 ndims is not 1'
       stop 1
    end if
    if (shape_in(1) /= nx*writerSize) then
       write(*,*) 'r64 shape_in read failed'
       stop 1
    end if
    deallocate(shape_in)

    call adios2_inquire_variable(variables(7), ioRead, "r64_2d", ierr)
    if (variables(7)%name /= 'r64_2d') then
       write(*,*) 'r64_2d not recognized'
       stop 1
    end if
    if (variables(7)%type /= adios2_type_dp) then
       write(*,*) 'r64_2d type not recognized'
       stop 1
    end if
    call adios2_variable_shape( shape_in, ndims, variables(7), ierr)
    if (ndims /= 2) then
       write(*,*) 'r64_2d ndims is not 2'
       stop 1
    end if
    if (shape_in(1) /= 2) then
       write(*,*) 'r64_2d shape_in(1) read failed'
       stop 1
    end if
    if (shape_in(2) /= nx*writerSize) then
       write(*,*) 'r64_2d shape_in(2) read failed'
       stop 1
    end if
    deallocate(shape_in)

    call adios2_inquire_variable(variables(8), ioRead, "r64_2d_rev", ierr)
    if (variables(8)%name /= 'r64_2d_rev') then
       write(*,*) 'r64_2d_rev not recognized'
       stop 1
    end if
    if (variables(8)%type /= adios2_type_dp) then
       write(*,*) 'r64_2d_rev type not recognized'
       stop 1
    end if
    call adios2_variable_shape(shape_in, ndims, variables(8), ierr)
    if (ndims /= 2) then
       write(*,*) 'r64_2d_rev ndims is not 2'
       stop 1
    end if
    if (shape_in(1) /= nx*writerSize) then
       write(*,*) 'r64_2d_rev shape_in(2) read failed'
       stop 1
    end if
    if (shape_in(2) /= 2) then
       write(*,*) 'r64_2d_rev shape_in(1) read failed'
       stop 1
    end if
    deallocate(shape_in)

    call adios2_inquire_variable(variables(10), ioRead, "c32", ierr)
    if (variables(10)%name /= 'c32') then
       write(*,*) 'c32 not recognized'
       stop 1
    end if
    if (variables(10)%type /= adios2_type_complex) then
       write(*,*) 'c32 type not recognized'
       stop 1
    end if
    call adios2_variable_shape(shape_in, ndims, variables(10), ierr)
    if (ndims /= 1) then
       write(*,*) 'c32 ndims is not 1'
       stop 1
    end if
    if (shape_in(1) /= nx*writerSize) then
       write(*,*) 'c32 shape_in read failed'
       stop 1
    end if
    deallocate(shape_in)

    call adios2_inquire_variable(variables(11), ioRead, "c64", ierr)
    if (variables(11)%name /= 'c64') then
       write(*,*) 'c64 not recognized'
       stop 1
    end if
    if (variables(11)%type /= adios2_type_complex_dp) then
       write(*,*) 'c64 type not recognized'
       stop 1
    end if
    call adios2_variable_shape(shape_in, ndims, variables(11), ierr)
    if (ndims /= 1) then
       write(*,*) 'c64 ndims is not 1'
       stop 1
    end if
    if (shape_in(1) /= nx*writerSize) then
       write(*,*) 'c64 shape_in read failed'
       stop 1
    end if
    deallocate(shape_in)

    call adios2_inquire_variable(variables(12), ioRead, "scalar_r64", ierr)
    if (variables(12)%name /= 'scalar_r64') then
       write(*,*) 'scalar_r64 not recognized'
       stop 1
    end if
    if (variables(12)%type /= adios2_type_dp) then
       write(*,*) 'scalar_r64 type not recognized'
       stop 1
    end if

    myStart = (writerSize * Nx / isize) * irank
    myLength = ((writerSize * Nx + isize - 1) / isize)

    if (myStart + myLength > writerSize * Nx) then
        myLength = writerSize * Nx - myStart
    end if
    allocate (in_I8(myLength));
    allocate (in_I16(myLength));
    allocate (in_I32(myLength));
    allocate (in_I64(myLength));
    allocate (in_R32(myLength));
    allocate (in_R64(myLength));
    allocate (in_C32(myLength));
    allocate (in_C64(myLength));
    allocate (in_R64_2d(2, myLength));
    allocate (in_R64_2d_rev(myLength, 2));

    start_dims(1) = myStart
    count_dims(1) = myLength

    start_dims2(1) = 0
    count_dims2(1) = 2
    start_dims2(2) = myStart
    count_dims2(2) = myLength

    start_dims3(1) = myStart
    count_dims3(1) = myLength
    start_dims3(2) = 0
    count_dims3(2) = 2

    call adios2_set_selection(variables(1), 1, start_dims, count_dims, ierr)
    call adios2_set_selection(variables(2), 1, start_dims, count_dims, ierr)
    call adios2_set_selection(variables(3), 1, start_dims, count_dims, ierr)
    call adios2_set_selection(variables(4), 1, start_dims, count_dims, ierr)
    call adios2_set_selection(variables(5), 1, start_dims, count_dims, ierr)
    call adios2_set_selection(variables(6), 1, start_dims, count_dims, ierr)
    call adios2_set_selection(variables(7), 2, start_dims2, count_dims2, ierr)
    call adios2_set_selection(variables(8), 2, start_dims3, count_dims3, ierr)
    call adios2_set_selection(variables(10), 1, start_dims, count_dims, ierr)
    call adios2_set_selection(variables(11), 1, start_dims, count_dims, ierr)

    call adios2_get(sstReader, variables(1), in_I8, ierr)
    call adios2_get(sstReader, variables(2), in_I16, ierr)
    call adios2_get(sstReader, variables(3), in_I32, ierr)
    call adios2_get(sstReader, variables(4), in_I64, ierr)
    call adios2_get(sstReader, variables(5), in_R32, ierr)
    call adios2_get(sstReader, variables(6), in_R64, ierr)
    call adios2_get(sstReader, variables(7), in_R64_2d, ierr)
    call adios2_get(sstReader, variables(8), in_R64_2d_rev, ierr)
    call adios2_get(sstReader, variables(10), in_C32, ierr)
    call adios2_get(sstReader, variables(11), in_C64, ierr)
    call adios2_get(sstReader, variables(12), in_scalar_R64, ierr)

    call adios2_end_step(sstReader, ierr)

    if(insteps == 0) then
      call ValidateTestData(myStart, myLength, 0)
    endif
    deallocate (in_I8);
    deallocate (in_I16);
    deallocate (in_I32);
    deallocate (in_I64);
    deallocate (in_R32);
    deallocate (in_R64);
    deallocate (in_C32);
    deallocate (in_C64);
    deallocate (in_R64_2d);
    deallocate (in_R64_2d_rev);
    insteps = insteps + 1
  end do

  call adios2_close(sstReader, ierr)

  ! Deallocates adios and calls its destructor
  call adios2_finalize(adios, ierr)
  if( adios%valid .eqv. .true. ) then
     write(*,*) 'Invalid adios2_finalize'
     stop 1
  end if


#if ADIOS2_USE_MPI
#ifdef CRAY_MPICH_VERSION
  call MPI_Barrier(MPI_COMM_WORLD)
#else
  call MPI_Finalize(ierr)
#endif
#endif


 end program TestSstRead
