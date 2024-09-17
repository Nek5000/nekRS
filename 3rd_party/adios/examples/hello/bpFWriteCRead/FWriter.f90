program FWriter
    use mpi
    use adios2
    implicit none

    integer(kind=8), dimension(2) :: shape_dims, start_dims, count_dims
    real, dimension(:,:), allocatable :: data
    integer :: i, j, inx, iny, irank, isize, ierr

    ! adios2 handlers
    type(adios2_adios):: adios
    type(adios2_io):: io
    type(adios2_variable):: var
    type(adios2_engine):: engine

    ! Launch MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

    ! Application variables
    inx = 3
    iny = 4
    allocate( data(inx, iny) )

    do j=1,iny
      do i=1,inx
        data(i,j) = irank * inx * iny +  (j-1)*inx + (i-1)
      end do
    end do

    ! Variable dimensions
    shape_dims(1) = isize * inx
    shape_dims(2) = iny

    start_dims(1) = irank * inx
    start_dims(2) = 0

    count_dims(1) = inx
    count_dims(2) = iny

    ! Create adios handler passing the communicator and error flag
    call adios2_init(adios, MPI_COMM_WORLD, ierr)

    ! Declare an IO process configuration inside adios
    call adios2_declare_io(io, adios, "FWriter", ierr)

    ! Defines a variable to be written in bp format
    call adios2_define_variable(var, io, "data2D", adios2_type_real, 2, &
                                shape_dims, start_dims, count_dims, &
                                adios2_constant_dims, ierr)

    ! Open in write mode, this launches an engine
    call adios2_open(engine, io, "FWriter.bp", adios2_mode_write, ierr)

    ! Put data contents to bp buffer
    call adios2_put(engine, var, data, ierr)

    ! Closes engine1 and deallocates it, becomes unreachable
    call adios2_close(engine, ierr)

    call adios2_finalize(adios, ierr)

    deallocate(data)

    call MPI_Finalize(ierr)

end program FWriter
