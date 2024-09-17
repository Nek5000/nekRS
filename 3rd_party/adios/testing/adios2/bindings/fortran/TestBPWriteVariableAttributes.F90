program TestBPWriteVariableAttributes
    use small_test_data
    use mpi
    use adios2
    implicit none


    type(adios2_adios) :: adios
    type(adios2_io) :: ioWrite
    type(adios2_engine) :: bpWriter
    type(adios2_variable) :: var
    type(adios2_attribute), dimension(14) :: attributes
!    type(adios2_attribute) :: failed_att

    integer :: ierr, i

    ! Launch MPI
    INTEGER provided

    ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)

    ! Create adios handler passing the communicator and error flag
    call adios2_init(adios, MPI_COMM_WORLD, ierr)

    !!!!!!!!!!!!!!!!!!!!!!!!! WRITER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Declare an IO process configuration inside adios
    call adios2_declare_io(ioWrite, adios, "ioWrite", ierr)

    ! create a global variable
    call adios2_define_variable(var, ioWrite, "myVar", adios2_type_integer4, &
                                ierr)

    ! Failed with ci/circleci: suse-pgi-openmpi (exceptions trigger abort)
!    call adios2_define_attribute(failed_att, ioWrite, 'att_String', &
!                                 'ADIOS2 String attribute', 'myVar2', '/', ierr)
!    if(ierr == 0) then 'myVar2 does not exist, should not create attribute att_String'
!    if(failed_att%valid .eqv. .true.) then
!        stop 'failed attribute must not exist '
!    end if

    do i=1,14
       if( attributes(i)%valid .eqv. .true. ) then
          write(*,*) 'Invalid attribute default'
          stop 1
       end if
    end do

    ! single value
    call adios2_define_attribute(attributes(1), ioWrite, 'att_String', &
                                 'ADIOS2 String attribute', var%name, '/', ierr)

    call adios2_define_attribute(attributes(2), ioWrite, 'att_i8', &
                                 data_I8(1), var%name, '/', ierr)

    call adios2_define_attribute(attributes(3), ioWrite, 'att_i16', &
                                 data_I16(1), var%name, '/', ierr)

    call adios2_define_attribute(attributes(4), ioWrite, 'att_i32', &
                                 data_I32(1), var%name, '/', ierr)

    call adios2_define_attribute(attributes(5), ioWrite, 'att_i64', &
                                 data_I64(1), var%name, '/', ierr)

    call adios2_define_attribute(attributes(6), ioWrite, 'att_r32', &
                                 data_R32(1), var%name, '/', ierr)

    call adios2_define_attribute(attributes(7), ioWrite, 'att_r64', &
                                 data_R64(1), var%name, '/', ierr)

    ! arrays
    call adios2_define_attribute(attributes(8), ioWrite, 'att_Strings_array', &
                                 data_Strings, 3, var%name, '/', ierr)

    call adios2_define_attribute(attributes(9), ioWrite, 'att_i8_array', &
                                 data_I8, 3, var%name, '/', ierr)

    call adios2_define_attribute(attributes(10), ioWrite, 'att_i16_array', &
                                 data_I16, 3, var%name, '/', ierr)

    call adios2_define_attribute(attributes(11), ioWrite, 'att_i32_array', &
                                 data_I32, 3, var%name, '/', ierr)

    call adios2_define_attribute(attributes(12), ioWrite, 'att_i64_array', &
                                 data_I64, 3, var%name, '/', ierr)

    ! uses / by default
    call adios2_define_attribute(attributes(13), ioWrite, 'att_r32_array', &
                                 data_R32, 3, var%name, ierr)

    call adios2_define_attribute(attributes(14), ioWrite, 'att_r64_array', &
                                 data_R64, 3, var%name, ierr)

    do i=1,14
       if( attributes(i)%valid .eqv. .false. ) then
          write(*,*) 'Invalid adios2_define_attribute'
          stop 1
       end if
    end do

    call adios2_open(bpWriter, ioWrite, "fvarattr_types.bp", adios2_mode_write, &
                     ierr)

    call adios2_put(bpWriter, var, 10, ierr)
    call adios2_close(bpWriter, ierr)

    call adios2_finalize(adios, ierr)

    call MPI_Finalize(ierr)

end program TestBPWriteVariableAttributes
