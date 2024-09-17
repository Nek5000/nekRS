! Helper functions for all examples
module decomp
contains

! random integer from {minv, minv+1, ..., maxv}
! including minv and maxv
function get_random(minv, maxv) result(n) 
    implicit none
    integer, intent(in) :: minv, maxv
    real :: r
    integer :: n
    call random_number(r)
    n = minv + FLOOR((maxv+1-minv)*r)  
end function get_random
    
! gather the local sizes of arrays and sum them up
! so that each process knows the global shape
! and its own offset in the global space
subroutine gather_decomp_1d(mysize, myshape, myoffset)
    use mpivars
    implicit none
    integer*8, intent(in) :: mysize
    integer*8, intent(out) :: myshape, myoffset
    integer*8, dimension(:), allocatable :: sizes

    allocate(sizes(nproc))
    call MPI_Allgather( mysize, 1, MPI_LONG_LONG, &
                        sizes, 1, MPI_LONG_LONG, &
                        app_comm, ierr)
    myshape = sum(sizes)
    myoffset = sum(sizes(1:rank))
    deallocate(sizes)
end subroutine gather_decomp_1d

subroutine decompose_1d(globalsize, myoffset, mysize)
    use mpivars
    implicit none
    integer*8, intent(in) :: globalsize
    integer*8, intent(out) :: myoffset, mysize
    integer*8 :: rem

    mysize = globalsize/nproc
    rem = globalsize-(nproc*mysize)
    if (rank < rem) then
        mysize = mysize + 1
        myoffset = rank*mysize 
    else
        myoffset = rank*mysize + rem
    endif
end subroutine decompose_1d

end module decomp

