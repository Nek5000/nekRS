program adios2_global_array_1d_read
    use mpivars
    use adios2
    implicit none

    ! ADIOS2 variables
    type(adios2_adios) :: adios

    ! MPI then ADIOS2 initialization
    call init_mpi(234)
    call adios2_init(adios, app_comm, ierr)

    call reader()

    ! ADIOS2 then MPI finalization
    call adios2_finalize(adios, ierr)
    call finalize_mpi()

contains

subroutine reader
    use mpivars
    use decomp
    use adios2
    implicit none

    character(len=256), parameter :: streamname = "adios2-global-array-1d-f.bp"

    type(adios2_io) :: io
    type(adios2_engine) :: engine
    type(adios2_variable) :: var_g
    integer :: step, istatus

    ! Application variables
    ! g = 1D distributed array, global shape and per-process size is fixed

    real*4, dimension(:), allocatable :: g
    integer :: ndims
    integer*8, dimension(:), allocatable :: fixed_shape
    integer*8, dimension(1) :: fixed_start, fixed_count

    call adios2_declare_io (io, adios, 'input', ierr)
    call adios2_open(engine, io, streamname, adios2_mode_read, ierr)
    if (ierr .ne. 0) then
        print '(" Failed to open stream: ",a)', streamname
        print '(" open stream ierr=: ",i0)', ierr
        return
    endif

    ! Reading steps
    step = 0
    do
        call adios2_begin_step(engine, adios2_step_mode_read, 0.0, istatus, ierr)
        if (ierr /= 0) then
            print '(" Failure when trying to get next step: ",a)', streamname
            exit
        endif
        if (istatus == adios2_step_status_end_of_stream) then
            ! Stream has terminated, no more steps are available
            !print '(" Input stream has terminated: ",a)', streamname
            exit
        endif

        ! Variable pointer MUST be retrieved every step, the reference
        ! will go invalid after adios2_end_step
        call adios2_inquire_variable(var_g, io, "GlobalArray", ierr )

        ! Get variable dimensions and do decomposition in the first step
        ! These don't change for the stream in this example
        if (step == 0) then
            ! fixed_shape is allocated in the next call
            call adios2_variable_shape(fixed_shape, ndims, var_g, ierr)

            call decompose_1d(fixed_shape(1), fixed_start(1), fixed_count(1))
            allocate(g(fixed_count(1)))

            write (*,100) "Read plan rank=", rank, &
            " global shape = ", fixed_shape(1),  &
            " local count = ", fixed_count(1),  &
            " offset = ", fixed_start(1) 
100 format (a,i2,a,i4,a,i1,a,i4)
        endif

        call adios2_set_selection(var_g, 1, fixed_start, fixed_count, ierr)
        call adios2_get(engine, var_g, g, ierr)        
        call adios2_end_step(engine, ierr)

        ! g[] is now filled with data AFTER adios2_end_step/adios2_perform_gets
        ! or should call adios2_get(engine, var_g, g, adios2_mode_sync, ierr)
        ! to get it immediately in the get() call

        step = step + 1
       enddo

    ! Close the output
    call adios2_close(engine, ierr)

    deallocate(g)
    deallocate(fixed_shape)

end subroutine reader


end program adios2_global_array_1d_read
