 program TestBPWriteTypes
     use small_test_data
     use mpi
     use adios2
     implicit none

     integer(kind=8), dimension(1) :: shape_dims, start_dims, count_dims
     integer :: inx, irank, isize, ierr, i, step_status

     type(adios2_adios) :: adios
     type(adios2_io) :: ioWrite, ioRead
     type(adios2_variable), dimension(12) :: variables
     type(adios2_engine) :: bpWriter, bpReader

     ! read handlers
     integer :: ndims
     integer(kind=8), dimension(:), allocatable :: shape_in

     ! Launch MPI
     INTEGER provided

     ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
     call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
     call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

     ! Application variables
     inx = 10

     ! Variable dimensions
     shape_dims(1) = isize*inx
     start_dims(1) = irank*inx
     count_dims(1) = inx

     ! Create adios handler passing the communicator and error flag
     call adios2_init(adios, MPI_COMM_WORLD, ierr)

     !!!!!!!!!!!!!!!!!!!!!!!!! WRITER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Declare an IO process configuration inside adios
     call adios2_declare_io(ioWrite, adios, "ioWrite", ierr)

     ! Defines a variable to be written in bp format
     call adios2_define_variable(variables(1), ioWrite, "var_I8", &
                                 adios2_type_integer1, 1, &
                                 shape_dims, start_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(2), ioWrite, "var_I16", &
                                 adios2_type_integer2, 1, &
                                 shape_dims, start_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(3), ioWrite, "var_I32", &
                                 adios2_type_integer4, 1, &
                                 shape_dims, start_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(4), ioWrite, "var_I64", &
                                 adios2_type_integer8, 1, &
                                 shape_dims, start_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(5), ioWrite, "var_R32", &
                                 adios2_type_real, 1, &
                                 shape_dims, start_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(6), ioWrite, "var_R64", &
                                 adios2_type_dp, 1, &
                                 shape_dims, start_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     ! Global variables
     call adios2_define_variable(variables(7), ioWrite, "gvar_I8", &
                                 adios2_type_integer1,  ierr)

     call adios2_define_variable(variables(8), ioWrite, "gvar_I16", &
                                 adios2_type_integer2,  ierr)

     call adios2_define_variable(variables(9), ioWrite, "gvar_I32", &
                                 adios2_type_integer4,  ierr)

     call adios2_define_variable(variables(10), ioWrite, "gvar_I64", &
                                 adios2_type_integer8,  ierr)

     call adios2_define_variable(variables(11), ioWrite, "gvar_R32", &
                                 adios2_type_real,  ierr)

     call adios2_define_variable(variables(12), ioWrite, "gvar_R64", &
                                 adios2_type_dp,  ierr)

     ! Open myVector_f.bp in write mode, this launches an engine
     call adios2_open(bpWriter, ioWrite, "ftypes.bp", adios2_mode_write, ierr)

     do i = 1, 3
         call adios2_begin_step(bpWriter, adios2_step_mode_append, -1.0, &
                               step_status, ierr)
         ! Put array contents to bp buffer, based on var1 metadata
         if (irank == 0 .and. i == 1) then
           call adios2_put(bpWriter, "gvar_I8", data_I8(1), ierr)
           call adios2_put(bpWriter, "gvar_I16", data_I16(1), ierr)
           call adios2_put(bpWriter, "gvar_I32", data_I32(1), ierr)
           call adios2_put(bpWriter, "gvar_I64", data_I64(1), ierr)
           call adios2_put(bpWriter, "gvar_R32", data_R32(1), ierr)
           call adios2_put(bpWriter, "gvar_R64", data_R64(1), ierr)
         end if

         call adios2_put(bpWriter, "var_I8", data_I8, ierr)
         call adios2_put(bpWriter, "var_I16", data_I16, ierr)
         call adios2_put(bpWriter, "var_I32", data_I32, ierr)
         call adios2_put(bpWriter, "var_I64", data_I64, ierr)
         call adios2_put(bpWriter, "var_R32", data_R32, ierr)
         call adios2_put(bpWriter, "var_R64", data_R64, ierr)
         call adios2_end_step(bpWriter, ierr)
     end do

     ! Closes engine1 and deallocates it, becomes unreachable
     call adios2_close(bpWriter, ierr)

     !!!!!!!!!!!!!!!!!!!!!!!!! READER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Declare io reader
     call adios2_declare_io(ioRead, adios, "ioRead", ierr)
     ! Open bpReader engine
     call adios2_open(bpReader, ioRead, "ftypes.bp", adios2_mode_read, ierr)

     call adios2_begin_step(bpReader, adios2_step_mode_read, -1., &
                            step_status, ierr)

     call adios2_inquire_variable(variables(1), ioRead, "var_I8", ierr)
     if (variables(1)%name /= 'var_I8') then
        write(*,*) 'var_I8 not recognized'
        stop 1
     end if
     if (variables(1)%type /= adios2_type_integer1) then
        write(*,*) 'var_I8 type not recognized'
        stop 1
     end if
     call adios2_variable_shape(shape_in, ndims, variables(1), ierr)
     if (ndims /= 1) then
        write(*,*) 'var_I8 ndims is not 1'
        stop 1
     end if
     if (shape_in(1) /= isize*inx) then
        write(*,*) 'var_I8 shape_in read failed'
        stop 1
     end if

     call adios2_inquire_variable(variables(2), ioRead, "var_I16", ierr)
     if (variables(2)%name /= 'var_I16') then
        write(*,*) 'var_I16 not recognized'
        stop 1
     end if
     if (variables(2)%type /= adios2_type_integer2) then
        write(*,*) 'var_I16 type not recognized'
        stop 1
     end if
     call adios2_variable_shape( shape_in, ndims,variables(2),ierr)
     if (ndims /= 1) then
        write(*,*) 'var_I16 ndims is not 1'
        stop 1
     end if
     if (shape_in(1) /= isize*inx) then
        write(*,*) 'var_I16 shape_in read failed'
        stop 1
     end if

     call adios2_inquire_variable(variables(3), ioRead, "var_I32", ierr)
     if (variables(3)%name /= 'var_I32') then
        write(*,*) 'var_I32 not recognized'
        stop 1
     end if
     if (variables(3)%type /= adios2_type_integer4) then
        write(*,*) 'var_I32 type not recognized'
        stop 1
     end if
     call adios2_variable_shape( shape_in, ndims, variables(3),ierr)
     if (ndims /= 1) then
        write(*,*) 'var_I32 ndims is not 1'
        stop 1
     end if
     if (shape_in(1) /= isize*inx) then
        write(*,*) 'var_I32 shape_in read failed'
        stop 1
     end if

     call adios2_inquire_variable(variables(4), ioRead, "var_I64", ierr)
     if (variables(4)%name /= 'var_I64') then
        write(*,*) 'var_I64 not recognized'
        stop 1
     end if
     if (variables(4)%type /= adios2_type_integer8) then
        write(*,*) 'var_I64 type not recognized'
        stop 1
     end if
     call adios2_variable_shape(shape_in, ndims, variables(4),ierr)
     if (ndims /= 1) then
        write(*,*) 'var_I64 ndims is not 1'
        stop 1
     end if
     if (shape_in(1) /= isize*inx) then
        write(*,*) 'var_I64 shape_in read failed'
        stop 1
     end if

     call adios2_inquire_variable(variables(5), ioRead, "var_R32", ierr)
     if (variables(5)%name /= 'var_R32') then
        write(*,*) 'var_R32 not recognized'
        stop 1
     end if
     if (variables(5)%type /= adios2_type_real) then
        write(*,*) 'var_R32 type not recognized'
        stop 1
     end if
     call adios2_variable_shape( shape_in, ndims, variables(5), ierr)
     if (ndims /= 1) then
        write(*,*) 'var_R32 ndims is not 1'
        stop 1
     end if
     if (shape_in(1) /= isize*inx) then
        write(*,*) 'var_R32 shape_in read failed'
        stop 1
     end if

     call adios2_inquire_variable(variables(6), ioRead, "var_R64", ierr)
     if (variables(6)%name /= 'var_R64') then
        write(*,*) 'var_R64 not recognized'
        stop 1
     end if
     if (variables(6)%type /= adios2_type_dp) then
        write(*,*) 'var_R64 type not recognized'
        stop 1
     end if
     call adios2_variable_shape( shape_in, ndims,  variables(6), ierr)
     if (ndims /= 1) then
        write(*,*) 'var_R64 ndims is not 1'
        stop 1
     end if
     if (shape_in(1) /= isize*inx) then
        write(*,*) 'var_R64 shape_in read failed'
        stop 1
     end if

     call adios2_end_step(bpReader, ierr)

     call adios2_close(bpReader, ierr)

     ! Deallocates adios and calls its destructor
     call adios2_finalize(adios, ierr)

     call MPI_Finalize(ierr)

 end program TestBPWriteTypes
