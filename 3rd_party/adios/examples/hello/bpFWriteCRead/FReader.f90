program FReader
    use mpi
    use adios2
    implicit none

    integer(kind=8), dimension(2) :: sel_start, sel_count
    real, dimension(:,:), allocatable :: data
    integer(kind=8) :: i, j
    integer :: irank, isize, ierr

    ! adios2 handlers
    type(adios2_adios):: adios
    type(adios2_io):: io
    type(adios2_variable):: var
    type(adios2_engine):: engine

    ! Launch MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

    if( irank == 0 ) then
        ! Create adios handler passing the communicator and error flag
        call adios2_init(adios, MPI_COMM_SELF, ierr)

        ! Declare an IO process configuration inside adios
        call adios2_declare_io(io, adios, "FReader", ierr)

        ! Open in write mode, this launches an engine
        call adios2_open(engine, io, "CppWriter.bp", adios2_mode_read, ierr)

        call adios2_inquire_variable(var, io, 'data2D', ierr)

        if( ierr == adios2_found ) then

            sel_start = (/ 0, 2 /)
            sel_count = (/ 3, 2 /)
            allocate( data( sel_count(1), sel_count(2) ) )

            call adios2_set_selection( var, 2, sel_start, sel_count, ierr )

            call adios2_get(engine, var, data, adios2_mode_sync, ierr)

            write(*,'(A,2(I2,A),A,2(I2,A),A)') 'Selection  &
                      & [ start = (', (sel_start(i),',',i=1,2) , ') &
                      &  count =  (', (sel_count(i),',',i=1,2) , ') ]'

            do j=1,sel_count(2)
              do i=1,sel_count(1)
                write(6,'(F3.0,A)', advance="no") data(i,j), ' '
              end do
              write(*,*)
            end do

            if( allocated(data) ) deallocate(data)

        end if

        call adios2_close(engine, ierr)
        call adios2_finalize(adios, ierr)

    end if

    call MPI_Finalize(ierr)

end program FReader
