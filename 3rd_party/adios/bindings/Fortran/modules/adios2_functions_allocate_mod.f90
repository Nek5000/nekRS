module adios2_functions_allocate_mod
    implicit none

    interface adios2_allocate

        ! 1D Array
        module procedure adios2_allocate_real_1d
        module procedure adios2_allocate_dp_1d
        module procedure adios2_allocate_complex_1d
        module procedure adios2_allocate_complex_dp_1d
        module procedure adios2_allocate_integer1_1d
        module procedure adios2_allocate_integer2_1d
        module procedure adios2_allocate_integer4_1d
        module procedure adios2_allocate_integer8_1d

        ! 2D Array
        module procedure adios2_allocate_real_2d
        module procedure adios2_allocate_dp_2d
        module procedure adios2_allocate_complex_2d
        module procedure adios2_allocate_complex_dp_2d
        module procedure adios2_allocate_integer1_2d
        module procedure adios2_allocate_integer2_2d
        module procedure adios2_allocate_integer4_2d
        module procedure adios2_allocate_integer8_2d

        ! 3D Array
        module procedure adios2_allocate_real_3d
        module procedure adios2_allocate_dp_3d
        module procedure adios2_allocate_complex_3d
        module procedure adios2_allocate_complex_dp_3d
        module procedure adios2_allocate_integer1_3d
        module procedure adios2_allocate_integer2_3d
        module procedure adios2_allocate_integer4_3d
        module procedure adios2_allocate_integer8_3d

        ! 4D Array
        module procedure adios2_allocate_real_4d
        module procedure adios2_allocate_dp_4d
        module procedure adios2_allocate_complex_4d
        module procedure adios2_allocate_complex_dp_4d
        module procedure adios2_allocate_integer1_4d
        module procedure adios2_allocate_integer2_4d
        module procedure adios2_allocate_integer4_4d
        module procedure adios2_allocate_integer8_4d

        ! 5D Array
        module procedure adios2_allocate_real_5d
        module procedure adios2_allocate_dp_5d
        module procedure adios2_allocate_complex_5d
        module procedure adios2_allocate_complex_dp_5d
        module procedure adios2_allocate_integer1_5d
        module procedure adios2_allocate_integer2_5d
        module procedure adios2_allocate_integer4_5d
        module procedure adios2_allocate_integer8_5d

        ! 6D Array
        module procedure adios2_allocate_real_6d
        module procedure adios2_allocate_dp_6d
        module procedure adios2_allocate_complex_6d
        module procedure adios2_allocate_complex_dp_6d
        module procedure adios2_allocate_integer1_6d
        module procedure adios2_allocate_integer2_6d
        module procedure adios2_allocate_integer4_6d
        module procedure adios2_allocate_integer8_6d

    end interface

contains

    ! 1D arrays
    subroutine adios2_allocate_real_1d(array, shp, ierr)
        real, dimension(:), allocatable, intent(out):: array
        integer(kind=8), dimension(1), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_dp_1d(array, shp, ierr)
        real(kind=8), dimension(:), allocatable, intent(out):: array
        integer(kind=8), dimension(1), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_1d(array, shp, ierr)
        complex, dimension(:), allocatable, intent(out):: array
        integer(kind=8), dimension(1), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_dp_1d(array, shp, ierr)
        complex(kind=8), dimension(:), allocatable, intent(out):: array
        integer(kind=8), dimension(1), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer1_1d(array, shp, ierr)
        integer(kind=1), dimension(:), allocatable, intent(out):: array
        integer(kind=8), dimension(1), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer2_1d(array, shp, ierr)
        integer(kind=2), dimension(:), allocatable, intent(out):: array
        integer(kind=8), dimension(1), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer4_1d(array, shp, ierr)
        integer(kind=4), dimension(:), allocatable, intent(out):: array
        integer(kind=8), dimension(1), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer8_1d(array, shp, ierr)
        integer(kind=8), dimension(:), allocatable, intent(out):: array
        integer(kind=8), dimension(1), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1)), stat=ierr)

    end subroutine

    ! 2D arrays
    subroutine adios2_allocate_real_2d(array, shp, ierr)
        real, dimension(:, :), allocatable, intent(out):: array
        integer(kind=8), dimension(2), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_dp_2d(array, shp, ierr)
        real(kind=8), dimension(:, :), allocatable, intent(out):: array
        integer(kind=8), dimension(2), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_2d(array, shp, ierr)
        complex, dimension(:, :), allocatable, intent(out):: array
        integer(kind=8), dimension(2), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_dp_2d(array, shp, ierr)
        complex(kind=8), dimension(:, :), allocatable, intent(out):: array
        integer(kind=8), dimension(2), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer1_2d(array, shp, ierr)
        integer(kind=1), dimension(:, :), allocatable, intent(out):: array
        integer(kind=8), dimension(2), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer2_2d(array, shp, ierr)
        integer(kind=2), dimension(:, :), allocatable, intent(out):: array
        integer(kind=8), dimension(2), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer4_2d(array, shp, ierr)
        integer(kind=4), dimension(:, :), allocatable, intent(out):: array
        integer(kind=8), dimension(2), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer8_2d(array, shp, ierr)
        integer(kind=8), dimension(:, :), allocatable, intent(out):: array
        integer(kind=8), dimension(2), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2)), stat=ierr)

    end subroutine

    ! 3D arrays
    subroutine adios2_allocate_real_3d(array, shp, ierr)
        real, dimension(:, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(3), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_dp_3d(array, shp, ierr)
        real(kind=8), dimension(:, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(3), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_3d(array, shp, ierr)
        complex, dimension(:, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(3), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_dp_3d(array, shp, ierr)
        complex(kind=8), dimension(:, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(3), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer1_3d(array, shp, ierr)
        integer(kind=1), dimension(:, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(3), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer2_3d(array, shp, ierr)
        integer(kind=2), dimension(:, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(3), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer4_3d(array, shp, ierr)
        integer(kind=4), dimension(:, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(3), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer8_3d(array, shp, ierr)
        integer(kind=8), dimension(:, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(3), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3)), stat=ierr)

    end subroutine

    ! 4D arrays
    subroutine adios2_allocate_real_4d(array, shp, ierr)
        real, dimension(:, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(4), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_dp_4d(array, shp, ierr)
        real(kind=8), dimension(:, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(4), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_4d(array, shp, ierr)
        complex, dimension(:, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(4), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_dp_4d(array, shp, ierr)
        complex(kind=8), dimension(:, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(4), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer1_4d(array, shp, ierr)
        integer(kind=1), dimension(:, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(4), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer2_4d(array, shp, ierr)
        integer(kind=2), dimension(:, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(4), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer4_4d(array, shp, ierr)
        integer(kind=4), dimension(:, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(4), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer8_4d(array, shp, ierr)
        integer(kind=8), dimension(:, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(4), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4)), stat=ierr)

    end subroutine

    ! 5D arrays
    subroutine adios2_allocate_real_5d(array, shp, ierr)
        real, dimension(:, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(5), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_dp_5d(array, shp, ierr)
        real(kind=8), dimension(:, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(5), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_5d(array, shp, ierr)
        complex, dimension(:, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(5), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_dp_5d(array, shp, ierr)
        complex(kind=8), dimension(:, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(5), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer1_5d(array, shp, ierr)
        integer(kind=1), dimension(:, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(5), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer2_5d(array, shp, ierr)
        integer(kind=2), dimension(:, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(5), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer4_5d(array, shp, ierr)
        integer(kind=4), dimension(:, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(5), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer8_5d(array, shp, ierr)
        integer(kind=8), dimension(:, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(5), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5)), stat=ierr)

    end subroutine

    ! 6D arrays
    subroutine adios2_allocate_real_6d(array, shp, ierr)
        real, dimension(:, :, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(6), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5), shp(6)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_dp_6d(array, shp, ierr)
        real(kind=8), dimension(:, :, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(6), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5), shp(6)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_6d(array, shp, ierr)
        complex, dimension(:, :, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(6), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5), shp(6)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_complex_dp_6d(array, shp, ierr)
        complex(kind=8), dimension(:, :, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(6), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5), shp(6)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer1_6d(array, shp, ierr)
        integer(kind=1), dimension(:, :, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(6), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5), shp(6)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer2_6d(array, shp, ierr)
        integer(kind=2), dimension(:, :, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(6), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5), shp(6)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer4_6d(array, shp, ierr)
        integer(kind=4), dimension(:, :, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(6), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5), shp(6)), stat=ierr)

    end subroutine

    subroutine adios2_allocate_integer8_6d(array, shp, ierr)
        integer(kind=8), dimension(:, :, :, :, :, :), allocatable, intent(out):: array
        integer(kind=8), dimension(6), intent(in):: shp
        integer, intent(out):: ierr

        if (allocated(array)) deallocate (array)
        allocate (array(shp(1), shp(2), shp(3), shp(4), shp(5), shp(6)), stat=ierr)

    end subroutine

end module
