program adios2_global_array_1d_write
    use mpivars
    use adios2
    implicit none
    integer, parameter :: numsteps = 5

    ! ADIOS2 variables
    type(adios2_adios) :: adios

    ! MPI then ADIOS2 initialization
    call init_mpi(123)
    call adios2_init(adios, app_comm, ierr)

    call writer()

    ! ADIOS2 then MPI finalization
    call adios2_finalize(adios, ierr)
    call finalize_mpi()

contains

subroutine writer
    use mpivars
    use decomp
    use adios2
    implicit none

    type(adios2_io) :: io
    type(adios2_engine) :: engine
    type(adios2_variable) :: var_g
    type(adios2_attribute) :: attr
    integer :: step

    ! Application variables
    ! g = 1D distributed array, 
    !   global shape and per-process size is fixed

    real*4, dimension(:), allocatable :: g
    character(80), parameter :: ga = "Global Array with fixed shape and decomposition"
    
    integer, parameter :: mincount = 2, maxcount = 5
    integer*8, dimension(1) :: fixed_shape, fixed_start, fixed_count

    fixed_count(1) = get_random(mincount, maxcount)
    allocate(g(fixed_count(1)))
    call gather_decomp_1d(fixed_count(1), fixed_shape(1), fixed_start(1))

    call adios2_declare_io (io, adios, 'output', ierr)

    call adios2_define_variable(var_g, io, "GlobalArray", &
                                adios2_type_real4, 1, &
                                fixed_shape, fixed_start, fixed_count, &
                                adios2_constant_dims, ierr)

    call adios2_define_attribute(attr, io, "GlobalArray/info", ga, ierr)

    call adios2_open(engine, io, "adios2-global-array-1d-f.bp", adios2_mode_write, ierr)

    write (*,100) "Decomp rank=", rank, &
    " global shape = ", fixed_shape(1),  &
    " local count = ", fixed_count(1),  &
    " offset = ", fixed_start(1) 
100 format (a,i2,a,i4,a,i1,a,i4)

    ! Computation/output loop
    do step=0,numsteps-1
        g = rank + (step+1)/100.0
        ! Output all data
        call adios2_begin_step(engine, adios2_step_mode_append, ierr)
        call adios2_put(engine, var_g, g, ierr);
        call adios2_end_step(engine, ierr)
    enddo

    ! Close the output
    call adios2_close(engine, ierr)

    deallocate(g)

    if (rank == 0) then
        write (*,*) "Try the following: "
        write (*,'(a,a,i4)') &
            "  bpls -la adios2-global-array-1d-f.bp ", &
            "GlobalArray -d -n ", fixed_shape(1)
        write (*,'(a,a,i4)') &
            "  bpls -la adios2-global-array-1d-f.bp ", &
            "GlobalArray -d -t -n ", fixed_shape(1)
        write (*,'(a)') &
            "  mpirun -n 2 ./adios2_basics_globalArray1DRead_f "
    endif
end subroutine writer


end program adios2_global_array_1d_write
