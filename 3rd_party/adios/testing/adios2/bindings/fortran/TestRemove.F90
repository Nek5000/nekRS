program TestRemove
    use small_test_data
#if ADIOS2_USE_MPI
    use mpi
#endif
    use adios2
    implicit none

    integer(kind=8), dimension(1) :: shape_dims, start_dims, count_dims
    integer :: inx, irank, isize, ierr
    logical :: res

    ! low-level
    type(adios2_adios) :: adios
    type(adios2_io) :: ioWrite
    type(adios2_variable), dimension(12) :: variables

#if ADIOS2_USE_MPI
    ! Launch MPI
    INTEGER provided

    ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)
#else
    irank = 0
    isize = 1
#endif

    ! Application variables
    inx = 10

    ! Variable dimensions
    shape_dims(1) = isize*inx
    start_dims(1) = irank*inx
    count_dims(1) = inx

    ! Create adios handler passing the communicator and error flag
#if ADIOS2_USE_MPI
    call adios2_init(adios, MPI_COMM_WORLD, ierr)
#else
    call adios2_init(adios, ierr)
#endif

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

     if (variables(1)%valid .eqv. .false. ) then
        write(*,*) 'var_I8 not defined'
        stop 1
     end if
     if (variables(2)%valid .eqv. .false. ) then
        write(*,*) 'var_I16 not defined'
        stop 1
     end if
     if (variables(3)%valid .eqv. .false. ) then
        write(*,*) 'var_I32 not defined'
        stop 1
     end if
     if (variables(4)%valid .eqv. .false. ) then
        write(*,*) 'var_I64 not defined'
        stop 1
     end if
     if (variables(5)%valid .eqv. .false. ) then
        write(*,*) 'var_R32 not defined'
        stop 1
     end if
     if (variables(6)%valid .eqv. .false. ) then
        write(*,*) 'var_R64 not defined'
        stop 1
     end if
     if (variables(7)%valid .eqv. .false. ) then
        write(*,*) 'gvar_I8 not defined'
        stop 1
     end if
     if (variables(8)%valid .eqv. .false. ) then
        write(*,*) 'gvar_I16 not defined'
        stop 1
     end if
     if (variables(9)%valid .eqv. .false. ) then
        write(*,*) 'gvar_I32 not defined'
        stop 1
     end if
     if (variables(10)%valid .eqv. .false. ) then
        write(*,*) 'gvar_I64 not defined'
        stop 1
     end if
     if (variables(11)%valid .eqv. .false. ) then
        write(*,*) 'gvar_R32 not defined'
        stop 1
     end if
     if (variables(12)%valid .eqv. .false. ) then
        write(*,*) 'gvar_IR64 not defined'
        stop 1
     end if

    ! remove piece
    call adios2_remove_variable(res, ioWrite, "gvar_R64", ierr)
    if( res .eqv. .false. ) then
       write(*,*) 'adios2_remove_variable failed'
       stop 1
    end if

    call adios2_inquire_variable(variables(12), ioWrite, "gvar_R64", ierr)
    if (variables(12)%valid .eqv. .true. ) then
       write(*,*) 'gvar_R64 found with inquire, not removed'
       stop 1
    end if

    ! remove all
    call adios2_remove_all_variables(ioWrite, ierr)

    call adios2_inquire_variable(variables(1), ioWrite, "var_I8", ierr)
    if (variables(1)%valid .eqv. .true. ) then
       write(*,*) 'var_I8 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(2), ioWrite, "var_I16", ierr)
    if (variables(2)%valid .eqv. .true.) then
       write(*,*) 'var_I16 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(3), ioWrite, "var_I32", ierr)
    if (variables(3)%valid .eqv. .true.) then
       write(*,*) 'var_I32 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(4), ioWrite, "var_I64", ierr)
    if (variables(4)%valid .eqv. .true.) then
       write(*,*) 'var_I64 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(5), ioWrite, "var_R32", ierr)
    if (variables(5)%valid .eqv. .true.) then
       write(*,*) 'var_R32 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(6), ioWrite, "var_R64", ierr)
    if (variables(6)%valid .eqv. .true.) then
       write(*,*) 'var_R64 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(7), ioWrite, "gvar_I8", ierr)
    if (variables(7)%valid .eqv. .true.) then
       write(*,*) 'gvar_I8 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(8), ioWrite, "gvar_I16", ierr)
    if (variables(8)%valid .eqv. .true.) then
       write(*,*) 'gvar_I16 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(9), ioWrite, "gvar_I32", ierr)
    if (variables(9)%valid .eqv. .true.) then
       write(*,*) 'gvar_I32 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(10), ioWrite, "gvar_I64", ierr)
    if (variables(10)%valid .eqv. .true.) then
       write(*,*) 'gvar_I64 found'
       stop 1
    end if

    call adios2_inquire_variable(variables(11), ioWrite, "gvar_R32", ierr)
    if (variables(11)%valid .eqv. .true.) then
       write(*,*) 'gvar_R32 found'
       stop 1
    end if

    call adios2_remove_io(res, adios, 'ioWrite', ierr)
    if( res .neqv. .true. ) then
       write(*,*) 'could not remove ioWrite'
       stop 1
    end if

    call adios2_at_io(ioWrite, adios, 'ioWrite', ierr)
    if( ioWrite%valid .eqv. .true. ) then
       write(*,*) 'did not remove ioWrite correctly'
       stop 1
    end if

    call adios2_finalize(adios, ierr)

#if ADIOS2_USE_MPI
    call MPI_Finalize(ierr)
#endif

end program TestRemove
