program adios2_values
    use mpivars
    use adios2
    implicit none
    integer, parameter :: numsteps = 5
    ! ADIOS2 variables
    type(adios2_adios) :: adios

    ! MPI then ADIOS2 initialization
    call init_mpi(12345)
    call adios2_init(adios, app_comm, ierr)

    call writer()
    call reader()

    ! ADIOS2 then MPI finalization
    call adios2_finalize(adios, ierr)
    call finalize_mpi()

contains

subroutine writer
    use mpivars
    use adios2
    implicit none

    type(adios2_variable) :: var_gc, var_gv, var_lc, var_lv, var_gs
    type(adios2_io) :: io
    type(adios2_engine) :: engine
    integer :: step

    ! Application variables
    ! gc = Global Constant (a single value from the entire application, once per run (e.g. NPROC))
    ! gv = Global Value (a single value from the entire application, changes over time)
    ! lc = Local Constant (one value per process, once per run (e.g. RANK))
    ! lv = Local Value (one value per process, changes over time)
    ! gs = a string, same as a global value
    integer :: gc, gv, lc, lv
    character(len=80) :: gs

    call adios2_declare_io (io, adios, 'Values', ierr)

    call adios2_define_variable(var_gc, io, "GlobalConstant", adios2_type_integer4, ierr)
    call adios2_define_variable(var_gv, io, "GlobalValue", adios2_type_integer4, ierr)
    ! Local values will show up in reading as a 1D array of nproc elements
    ! the write side definition is quite cumbersome in Fortran :-(
    ! We have to define it almost like a distributed global array with a special
    ! dimension value to indicate its type
    call adios2_define_variable(var_lc, io, "LocalConstant", adios2_type_integer4, &
                                1, (/ adios2_local_value_dim /), &
                                adios2_null_dims, &
                                adios2_null_dims, &
                                adios2_constant_dims, ierr)

    call adios2_define_variable(var_lv, io, "LocalValue", adios2_type_integer4, &
                                1, (/ adios2_local_value_dim /), &
                                adios2_null_dims, &
                                adios2_null_dims, &
                                adios2_constant_dims, ierr)
                                
    call adios2_define_variable(var_gs, io, "GlobalString", adios2_type_string, ierr)

    call adios2_open(engine, io, "adios2-values-f.bp", adios2_mode_write, ierr)

    ! Computation/output loop
    gc = nproc
    lc = rank
    do step=0,numsteps-1
        gv = step
        lv = nproc*(step)+rank
        write (gs,'(a,i3)') "This is step ", step

        call adios2_begin_step(engine, adios2_step_mode_append, ierr)
        if (step == 0) then
            call adios2_put(engine, var_lc, lc, ierr);
            if (rank == 0) then
                ! could be written from every process but it is useless
                call adios2_put(engine, var_gc, gc, ierr);
            endif
        endif
        if (rank == 0) then
            ! could be written from every process but it is useless
            call adios2_put(engine, var_gv, gv, ierr);
            call adios2_put(engine, var_gs, gs, ierr);
        endif
        call adios2_put(engine, var_lv, lv, ierr);
        call adios2_end_step(engine, ierr)
    enddo

    ! Close the output
    call adios2_close(engine, ierr)

end subroutine writer


subroutine reader
    use mpivars
    implicit none

    type(adios2_variable) :: var_gc, var_gv, var_lc, var_lv
    type(adios2_io) :: io
    type(adios2_engine) :: engine
    integer*8 :: numsteps, i
    integer*4 :: gc
    integer*4, dimension(:), allocatable :: gvs, lcs
    integer*4, dimension(:,:), allocatable :: lvs
    character(len=80)::fmt
    integer*8, dimension(:), allocatable :: shape_lv, shape_lc
    integer :: ndims_lv, ndims_lc

    ! Note, every process reads everything in this example
    
    call adios2_declare_io(io, adios, "ValuesInput", ierr)
    call adios2_open(engine, io, "adios2-values-f.bp", adios2_mode_readRandomAccess, MPI_COMM_SELF, ierr)
    
    call adios2_inquire_variable(var_gc, io, "GlobalConstant", ierr)
    call adios2_get(engine, var_gc, gc , ierr)

    call adios2_inquire_variable(var_gv, io, "GlobalValue", ierr)
    call adios2_variable_steps(numsteps, var_gv, ierr)
    call adios2_set_step_selection(var_gv, 0_8, numsteps, ierr)
    allocate(gvs(numsteps))
    call adios2_get(engine, var_gv, gvs , ierr)

    ! Read Local Values and Local Constants as a 1D array
    ! shape array is allocated inside adios2_variable_shape()

    call adios2_inquire_variable(var_lc, io, "LocalConstant", ierr)
    call adios2_variable_shape(shape_lc, ndims_lc, var_lc, ierr)
    allocate(lcs(shape_lc(1)))
    call adios2_get(engine, var_lc, lcs , ierr)

    call adios2_inquire_variable(var_lv, io, "LocalValue", ierr)
    call adios2_variable_shape(shape_lv, ndims_lv, var_lv, ierr)
    call adios2_set_step_selection(var_lv, 0_8, numsteps, ierr)
    allocate(lvs(shape_lv(1),numsteps))
    call adios2_get(engine, var_lv, lvs , ierr)

    call adios2_close(engine, ierr)

    ! By default, get()s are deferred and content is available AFTER
    ! adios2_close() or adios2_perform_gets()
    ! Use adios2_mode_sync option in adios2_get() to get the content immediately

    if (rank == 0) then
        write(*,'("Number of steps in file = ",i5)') numsteps
        write(*,'("GlobalConstant = ", i5)') gc

        write(fmt,'(a,i5,a)') '(a18,',numsteps,'i4,a2)'
        !write(*,'(a)') fmt
        write(*,fmt) "GlobalValue(s) = [", gvs, " ]"

        write(fmt,'(a,i5,a)') '(a20,',shape_lc(1),'i4,a2)'
        !write(*,'(a)') fmt
        write(*,fmt) "LocalConstant(s) = [", lcs, " ]"

        write(fmt,'(a,i5,a)') '(a6,i3,a4,',shape_lv(1),'i4)'
        !write(*,'(a)') fmt
        write(*,'(a)') "LocalValues = ["
        do i = 1, numsteps
            write(*,fmt) "  step", i-1, ":", lvs(:,i)
        enddo
        write(*,*) "             ]"
    endif

    deallocate(gvs, lcs, lvs)
end subroutine reader

end program adios2_values
