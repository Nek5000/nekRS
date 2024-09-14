! Distributed under the OSI-approved Apache License, Version 2.0.  See
! accompanying file Copyright.txt for details.
!
! SmallTestData_mod.f90 : small Fortran 90 arrays data for tests
!
!  Created on: Aug 9, 2017
!      Author: William F Godoy godoywf@ornl.gov
!

module sst_test_data
    implicit none

    integer, parameter :: Nx = 10

    integer(kind=1), dimension(10) :: data_I8
    integer(kind=2), dimension(10) :: data_I16
    integer(kind=4), dimension(10) :: data_I32
    integer(kind=8), dimension(10) :: data_I64
    real(kind=4), dimension(10) :: data_R32
    real(kind=8), dimension(10) :: data_R64
    complex(kind=4), dimension(10) :: data_C32
    complex(kind=8), dimension(10) :: data_C64
    real (kind=8), dimension(2, 10) :: data_R64_2d
    real (kind=8), dimension(10, 2) :: data_R64_2d_rev
    real (kind=8) :: data_scalar_R64

    integer(kind=1), dimension(:), allocatable :: in_I8
    integer(kind=2), dimension(:), allocatable :: in_I16
    integer(kind=4), dimension(:), allocatable :: in_I32
    integer(kind=8), dimension(:), allocatable :: in_I64
    real(kind=4), dimension(:), allocatable :: in_R32
    real(kind=8), dimension(:), allocatable :: in_R64
    complex(kind=4), dimension(:), allocatable :: in_C32
    complex(kind=8), dimension(:), allocatable :: in_C64
    real (kind=8), dimension(:,:), allocatable :: in_R64_2d
    real (kind=8), dimension(:,:), allocatable :: in_R64_2d_rev
    real (kind=8) :: in_scalar_R64

    contains
    subroutine GenerateTestData(step, rank)
      INTEGER, INTENT(IN) :: step, rank


      integer (kind=8) :: i, j
      j =  rank * Nx * 10 + step;
      data_scalar_r64 = (step + 1) * 1.5D0;
      do i = 1, Nx
         data_I8(i) = INT(j + 10 * (i-1), 1);
         data_I16(i) = INT(j + 10 * (i-1), 2);
         data_I32(i) = INT(j + 10 * (i-1), 4);
         data_I64(i) = INT(j + 10 * (i-1), 8);
         data_R32(i) = INT(j + 10 * (i-1), 4);
         data_R64(i) = INT(j + 10 * (i-1), 8);
         data_C32(i) = cmplx(-(j + 10 * (i-1)), (j + 10 * (i-1)));
         data_C64(i) = cmplx(-(j + 10 * (i-1)), (j + 10 * (i-1)));
         data_R64_2d(1,i) = (j + 10 * (i-1));
         data_R64_2d(2,i) = 10000 + (j + 10 * (i-1));
         data_R64_2d_rev(i,1) = (j + 10 * (i-1));
         data_R64_2d_rev(i,2) = 10000 + (j + 10 * (i-1));
      end do

    end subroutine GenerateTestData

    subroutine ValidateTestData(start, length, step)
      INTEGER, INTENT(IN) :: start, length, step

      integer (kind=8) :: i
      do i = 1, length
         if (in_scalar_R64 /= (step + 1) * 1.5D0) then
            stop 'scalar_r64 value failed'
         end if
         if (in_I8(i) /= INT(((i - 1 + start)* 10 + step), 1)) then
            stop 'data_I8 value failed'
         end if
         if (in_I16(i) /= (i - 1 + start)* 10 + step) then
            stop 'data_I16 value failed'
         end if
         if (in_I32(i) /= (i - 1 + start)* 10 + step) then
            stop 'data_I32 value failed'
         end if
         if (in_I64(i) /= (i - 1 + start)* 10 + step) then
            stop 'data_I64 value failed'
         end if
         if (in_R32(i) /= INT((i - 1 + start)* 10 + step, 4)) then
            stop 'data_R32 value failed'
         end if
         if (in_R64(i) /= (i - 1 + start)* 10 + step) then
            stop 'data_R64 value failed'
         end if
         if (in_R64_2d(1, i) /= (i - 1 + start)* 10 + step) then
            stop 'data_R64 value failed'
         end if
         if (in_R64_2d(2, i) /= 10000 + (i - 1 + start)* 10 + step) then
            stop 'data_R64 value failed'
         end if
         if (in_R64_2d_rev(i, 1) /= (i - 1 + start)* 10 + step) then
            stop 'data_R64_rev value failed'
         end if
         if (in_R64_2d_rev(i, 2) /= 10000 + (i - 1 + start)* 10 + step) then
            stop 'data_R64_rev value failed'
         end if
      end do

    end subroutine ValidateTestData

  end module sst_test_data
