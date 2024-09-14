!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
! accompanying file Copyright.txt for details.
!
! inlineMixedLang.f90: example borrowed from bpFWriteCRead example, but using
! the inline engine. inlineMixedLang.f90 creates data and uses inline writer
! for making that data available in C++ functions (inlineMixedLang.cpp) that
! reads the data.
!
! Created on: Oct 20, 2021
!     Author: Caitlin Ross caitlin.ross@kitware.com
!

program InlineMixedLang
    use adios2
    implicit none

    interface
        subroutine open_reader(io)
            import adios2_io
            type(adios2_io), intent(in) :: io
        end subroutine

        subroutine close_reader()
        end subroutine

        subroutine analyze_data()
        end subroutine
    end interface

    integer(kind=8), dimension(2) :: shape_dims, start_dims, count_dims
    real, dimension(:,:), allocatable :: data
    integer :: i, j, inx, iny, istep, ierr

    ! adios2 handlers
    type(adios2_adios):: adios
    type(adios2_io):: io
    type(adios2_variable):: var
    type(adios2_engine):: writer

    ! Application variables
    inx = 3
    iny = 4
    allocate( data(inx, iny)  )

    ! Variable dimensions
    shape_dims(1) = inx
    shape_dims(2) = iny

    start_dims(1) = 0
    start_dims(2) = 0

    count_dims(1) = inx
    count_dims(2) = iny

    ! Create adios handler
    call adios2_init(adios, ierr)

    ! Declare an IO process configuration inside adios
    call adios2_declare_io(io, adios, "Inline-mixed", ierr)

    ! Defines a variable to be "written" (with inline engine, pointer to data
    ! is passed to the reader)
    call adios2_define_variable(var, io, "data2D", adios2_type_real, 2, &
                                shape_dims, start_dims, count_dims, &
                                adios2_constant_dims, ierr)

    ! set the engine to be Inline
    call adios2_set_engine(io, "Inline", ierr)

    ! Open in write mode, this launches an engine
    call adios2_open(writer, io, "writer", adios2_mode_write, ierr)

    ! first call to a C++ function to open the inline reader
    call open_reader(io)

    do istep=1,5
      ! begin the step for writing
      call adios2_begin_step(writer, adios2_step_mode_append, ierr)

      ! generate data for this time step
      write(*,*) "fortran step ", istep
      do j=1,iny
        do i=1,inx
          data(i,j) = (inx * iny +  (j-1)*inx + (i-1)) * istep
        end do
        write(*,*) (data(i,j), i=1,inx)
      end do

      ! Put data contents and end the step for the writer.
      call adios2_put(writer, var, data, ierr)
      call adios2_end_step(writer, ierr)

      ! analyze data with C++ function
      call analyze_data()
    end do

    ! Closes writer and deallocates it, becomes unreachable
    call adios2_close(writer, ierr)

    ! Closes reader on C++ side
    call close_reader()

    call adios2_finalize(adios, ierr)

    deallocate(data)

end program InlineMixedLang
