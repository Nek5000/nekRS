!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_attribute_mod.f90 : ADIOS2 Fortran bindings for overloaded
!                             adios2_attribute_data subroutines
!   Created on: Dec 10, 2018
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_attribute_data_mod
    use adios2_attribute_mod
    implicit none

    interface adios2_attribute_data

        ! Single Value
        module procedure adios2_attribute_data_string
        module procedure adios2_attribute_data_real
        module procedure adios2_attribute_data_dp
        module procedure adios2_attribute_data_integer1
        module procedure adios2_attribute_data_integer2
        module procedure adios2_attribute_data_integer4
        module procedure adios2_attribute_data_integer8

        ! 1D Array
        module procedure adios2_attribute_data_string_1d
        module procedure adios2_attribute_data_real_1d
        module procedure adios2_attribute_data_dp_1d
        module procedure adios2_attribute_data_integer1_1d
        module procedure adios2_attribute_data_integer2_1d
        module procedure adios2_attribute_data_integer4_1d
        module procedure adios2_attribute_data_integer8_1d

    end interface

    external adios2_attribute_data_f2c
    external adios2_attribute_value_f2c
contains

    ! Single value
    subroutine adios2_attribute_data_string(data, attribute, ierr)
        character*(*), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr

        call adios2_attribute_check_type(attribute, adios2_type_string, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_value_f2c(data, attribute%f2c, ierr)
        end if

    end subroutine


    subroutine adios2_attribute_data_real(data, attribute, ierr)
        real, intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr

        call adios2_attribute_check_type(attribute, adios2_type_real, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_value_f2c(data, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_dp(data, attribute, ierr)
        real(kind=8), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr

        call adios2_attribute_check_type(attribute, adios2_type_dp, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_value_f2c(data, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_integer1(data, attribute, ierr)
        integer(kind=1), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr

        call adios2_attribute_check_type(attribute, adios2_type_integer1, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_value_f2c(data, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_integer2(data, attribute, ierr)
        integer(kind=2), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr

        call adios2_attribute_check_type(attribute, adios2_type_integer2, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_value_f2c(data, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_integer4(data, attribute, ierr)
        integer(kind=4), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr

        call adios2_attribute_check_type(attribute, adios2_type_integer4, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_value_f2c(data, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_integer8(data, attribute, ierr)
        integer(kind=8), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr

        call adios2_attribute_check_type(attribute, adios2_type_integer8, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_value_f2c(data, attribute%f2c, ierr)
        end if

    end subroutine

    ! 1D Array
    subroutine adios2_attribute_data_string_1D(data, attribute, ierr)
        character*(*), dimension(:), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr
        ! local
        integer:: i
        character(len=adios2_string_array_element_max_size), &
            dimension(attribute%length):: dataMax

        call adios2_attribute_check_type(attribute, adios2_type_string, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_data_f2c(dataMax, attribute%length, &
                                           attribute%f2c, ierr)
        end if

        do i=1, attribute%length
            data(i) = trim(dataMax(i))
        end do

    end subroutine

    subroutine adios2_attribute_data_real_1d(data, attribute, ierr)
        real, dimension(:), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr
        ! local
        integer :: length

        call adios2_attribute_check_type(attribute, adios2_type_real, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_data_f2c(data, length, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_dp_1d(data, attribute, ierr)
        real(kind=8), dimension(:), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr
        ! local
        integer :: length

        call adios2_attribute_check_type(attribute, adios2_type_dp, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_data_f2c(data, length, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_integer1_1d(data, attribute, ierr)
        integer(kind=1), dimension(:), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr
        ! local
        integer :: length

        call adios2_attribute_check_type(attribute, adios2_type_integer1, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_data_f2c(data, length, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_integer2_1d(data, attribute, ierr)
        integer(kind=2), dimension(:), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr
        ! local
        integer :: length

        call adios2_attribute_check_type(attribute, adios2_type_integer2, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_data_f2c(data, length, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_integer4_1d(data, attribute, ierr)
        integer(kind=4), dimension(:), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr
        ! local
        integer :: length

        call adios2_attribute_check_type(attribute, adios2_type_integer4, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_data_f2c(data, length, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_data_integer8_1d(data, attribute, ierr)
        integer(kind=8), dimension(:), intent(out):: data
        type(adios2_attribute), intent(in):: attribute
        integer, intent(out):: ierr
        ! local
        integer :: length

        call adios2_attribute_check_type(attribute, adios2_type_integer8, &
                                        'attribute_data', ierr)
        if (ierr == 0) then
            call adios2_attribute_data_f2c(data, length, attribute%f2c, ierr)
        end if

    end subroutine

end module
