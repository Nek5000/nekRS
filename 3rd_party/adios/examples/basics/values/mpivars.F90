module mpivars
    implicit none
    include 'mpif.h'

    integer:: app_comm, rank, nproc
    integer:: wrank, wnproc
    integer:: ierr

contains

subroutine init_mpi(color)
    integer, intent(in):: color
    call MPI_Init(ierr)
    ! World comm spans all applications started with the same mpirun command 
    call MPI_Comm_rank(MPI_COMM_WORLD, wrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, wnproc, ierr)

    ! Have to split and create a 'world' communicator for this app only
    ! color must be unique for each application 
    call MPI_Comm_split (MPI_COMM_WORLD, color, wrank, app_comm, ierr)
    call MPI_Comm_rank (app_comm, rank, ierr)
    call MPI_Comm_size (app_comm, nproc , ierr)
end subroutine init_mpi

subroutine finalize_mpi()
    call MPI_Finalize(ierr)
end subroutine finalize_mpi

end module mpivars

