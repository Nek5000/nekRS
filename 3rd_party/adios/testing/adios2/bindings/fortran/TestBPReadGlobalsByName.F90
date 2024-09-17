program TestBPReadGlobalsByName
    use mpi
    use adios2
    implicit none

    integer :: irank, isize, ierr, step_status
    integer(kind=4) :: diag_1d_nsp
    real(kind=8) :: sml_outpsi

    ! adios2 handlers
    type(adios2_adios):: adios
    type(adios2_io):: ioWrite, ioRead
    type(adios2_variable):: var
    type(adios2_engine):: writer, reader

    ! Launch MPI
    INTEGER provided

    ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

    ! Start adios2
    call adios2_init(adios, MPI_COMM_SELF, ierr)

    ! writer
    call adios2_declare_io(ioWrite, adios, "FWriter", ierr)

    call adios2_define_variable(var, ioWrite, "diag_1d_nsp", &
                                adios2_type_integer4, ierr)
    call adios2_define_variable(var, ioWrite, "sml_outpsi", adios2_type_dp, &
                                ierr)

    call adios2_open(writer, ioWrite, "f0analysis_static.bp", &
                     adios2_mode_write, ierr)
    call adios2_put(writer, "diag_1d_nsp", 1, ierr)
    call adios2_put(writer, "sml_outpsi", 0.295477_8, ierr)
    call adios2_close(writer, ierr)

    if(ierr /= 0) then
       write(*,*) 'Problems writing'
       stop 1
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! reader
    call adios2_declare_io(ioRead, adios, "FReader", ierr)
    call adios2_open(reader, ioRead, "f0analysis_static.bp", adios2_mode_read, ierr)
    call adios2_begin_step(reader, adios2_step_mode_read, -1., &
                           step_status, ierr)

    call adios2_get(reader, "diag_1d_nsp", diag_1d_nsp, ierr)
    call adios2_get(reader, "sml_outpsi", sml_outpsi, ierr)
    call adios2_end_step(reader, ierr)

    call adios2_close(reader, ierr)

    call adios2_finalize(adios, ierr)

    !! Expect 1 and 0.295477 for diag_1d_nsp and sml_outpsi
    print *, irank, 'diag_1d_nsp', diag_1d_nsp
    print *, irank, 'sml_outpsi', sml_outpsi

    if( diag_1d_nsp /= 1 ) then
       write(*,*) 'diag_1d_nsp is not 1'
       stop 1
    end if

    if( sml_outpsi /= 0.295477_8 ) then
       write(*,*) 'sml_outpsi is not 0.295477'
       stop 1
    end if

    call MPI_Finalize(ierr)

end program TestBPReadGlobalsByName
