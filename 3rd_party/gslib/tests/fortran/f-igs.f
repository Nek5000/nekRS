      program figs
      implicit none

      include 'mpif.h'

      integer npmax
      parameter(npmax=16)

      integer ierror,handle,hwait,np,me,i,neighbors,count
      integer*8 id(npmax)

      real*8 answer(npmax),u(npmax)

      call mpi_init(ierror)
      call mpi_comm_size(mpi_comm_world,np,ierror)
      call mpi_comm_rank(mpi_comm_world,me,ierror)

      count=1
      if(me.gt.0) then
        id(count)=me
        count=count+1
      endif
      id(count)=me+1
      count=count+1
      if(me.lt.(np-1)) then
        id(count)=me+2
        count=count+1
      endif

      neighbors=count-1
!     gs_pairwise
      call fgslib_gs_setup_pick(handle,id,neighbors,mpi_comm_world,np,1)

      if(np.eq.1) then
        answer(1)=1.0
      else
        answer(1)=2.0
        answer(np)=2.0
        do i=2,np-1
          answer(i)=3.0
        enddo
      endif

      do i=1,neighbors
        u(i)=1.0
      enddo

      call fgslib_igs_op(handle,u,1,1,0,hwait)
      call fgslib_gs_op_wait(hwait)

      do i=1,neighbors
        if(abs(u(i)-answer(id(i)))>1e-16) then
          write(6,*) 'igs_op test failed'
        endif
      enddo

      call mpi_finalize(ierror)

      end
