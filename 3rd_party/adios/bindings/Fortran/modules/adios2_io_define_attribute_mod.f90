!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_io_define_attribute_mod.f90 : ADIOS2 Fortran bindings for IO class
!  overloaded (C++ template) function adios2_define_attribute
!
!   Created on: Dec 15, 2017
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_io_define_attribute_mod
    use adios2_parameters_mod
    implicit none

    interface adios2_define_attribute

        ! Single value
        module procedure adios2_define_attribute_string
        module procedure adios2_define_attribute_real
        module procedure adios2_define_attribute_dp
        module procedure adios2_define_attribute_integer1
        module procedure adios2_define_attribute_integer2
        module procedure adios2_define_attribute_integer4
        module procedure adios2_define_attribute_integer8

        ! 1D Array
        module procedure adios2_define_attribute_string_1d
        module procedure adios2_define_attribute_real_1d
        module procedure adios2_define_attribute_dp_1d
        module procedure adios2_define_attribute_integer1_1d
        module procedure adios2_define_attribute_integer2_1d
        module procedure adios2_define_attribute_integer4_1d
        module procedure adios2_define_attribute_integer8_1d

        ! variable attribute
        ! Single value
        module procedure adios2_define_variable_attribute_string
        module procedure adios2_define_variable_attribute_real
        module procedure adios2_define_variable_attribute_dp
        module procedure adios2_define_variable_attribute_integer1
        module procedure adios2_define_variable_attribute_integer2
        module procedure adios2_define_variable_attribute_integer4
        module procedure adios2_define_variable_attribute_integer8

        ! 1D Array
        module procedure adios2_define_variable_attribute_string_1d
        module procedure adios2_define_variable_attribute_real_1d
        module procedure adios2_define_variable_attribute_dp_1d
        module procedure adios2_define_variable_attribute_integer1_1d
        module procedure adios2_define_variable_attribute_integer2_1d
        module procedure adios2_define_variable_attribute_integer4_1d
        module procedure adios2_define_variable_attribute_integer8_1d

        ! Single value
        module procedure adios2_def_var_att_sep_string
        module procedure adios2_def_var_att_sep_real
        module procedure adios2_def_var_att_sep_dp
        module procedure adios2_def_var_att_sep_integer1
        module procedure adios2_def_var_att_sep_integer2
        module procedure adios2_def_var_att_sep_integer4
        module procedure adios2_def_var_att_sep_integer8

        ! 1D Array
        module procedure adios2_def_var_att_sep_string_1d
        module procedure adios2_def_var_att_sep_real_1d
        module procedure adios2_def_var_att_sep_dp_1d
        module procedure adios2_def_var_att_sep_integer1_1d
        module procedure adios2_def_var_att_sep_integer2_1d
        module procedure adios2_def_var_att_sep_integer4_1d
        module procedure adios2_def_var_att_sep_integer8_1d

    end interface
    external adios2_define_vattr_f2c
    external adios2_define_vattr_array_f2c

contains

    subroutine adios2_define_attribute_string(attribute, io, name, value, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        character*(*), intent(in):: value
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_string, &
                                     TRIM(ADJUSTL(value))//char(0), char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = name
            attribute%type = adios2_type_string
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_attribute_real(attribute, io, name, value, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real, intent(in):: value
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_real, value, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = name
            attribute%type = adios2_type_real
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_attribute_dp(attribute, io, name, value, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(kind=8), intent(in):: value
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_dp, value, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = name
            attribute%type = adios2_type_dp
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_attribute_integer1(attribute, io, name, value, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=1), intent(in):: value
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer1, value, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = name
            attribute%type = adios2_type_integer1
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_attribute_integer2(attribute, io, name, value, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=2), intent(in):: value
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer2, value, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = name
            attribute%type = adios2_type_integer2
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_attribute_integer4(attribute, io, name, value, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=4), intent(in):: value
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer4, value, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = name
            attribute%type = adios2_type_integer4
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_attribute_integer8(attribute, io, name, value, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=8), intent(in):: value
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer8, value, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = name
            attribute%type = adios2_type_integer8
            attribute%length = 1
        end if
    end subroutine

    ! 1D
    subroutine adios2_define_attribute_string_1d(attribute, io, name, &
                                                 data, length, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        character*(*), dimension(:), intent(in):: data
        integer, intent(in) :: length
        integer, intent(out) :: ierr

        ! local data with zero terminated character
        character(len=adios2_string_array_element_max_size), &
            dimension(length):: data_null_terminated

        integer :: i

        do i = 1, length
            data_null_terminated(i) = TRIM(ADJUSTL(data(i)))//char(0)
        end do

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_string, &
                                           data_null_terminated, length, char(0), char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = name
            attribute%type = adios2_type_string
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_attribute_real_1d(attribute, io, name, data, &
                                               length, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real, dimension(:), intent(in):: data
        integer, intent(in) :: length
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_real, data, length, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = name
            attribute%type = adios2_type_real
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_attribute_dp_1d(attribute, io, name, data, &
                                             length, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(kind=8), dimension(:), intent(in):: data
        integer, intent(in) :: length
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_dp, data, length, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = name
            attribute%type = adios2_type_dp
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_attribute_integer1_1d(attribute, io, name, &
                                                   data, length, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=1), dimension(:), intent(in):: data
        integer, intent(in) :: length
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer1, &
                                           data, length, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = name
            attribute%type = adios2_type_integer1
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_attribute_integer2_1d(attribute, io, name, &
                                                   data, length, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=2), dimension(:), intent(in):: data
        integer, intent(in) :: length
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer2, &
                                           data, length, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = name
            attribute%type = adios2_type_integer2
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_attribute_integer4_1d(attribute, io, name, &
                                                   data, length, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=4), dimension(:), intent(in):: data
        integer, intent(in) :: length
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer4, data, length, char(0), char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = name
            attribute%type = adios2_type_integer4
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_attribute_integer8_1d(attribute, io, name, &
                                                   data, length, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=8), dimension(:), intent(in):: data
        integer, intent(in) :: length
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer8, &
                                           data, length, char(0), char(0), ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = name
            attribute%type = adios2_type_integer8
            attribute%length = length
        end if
    end subroutine

    ! Variable attribute
    subroutine adios2_define_variable_attribute_string(attribute, io, name, &
                                                       value, variable_name, &
                                                       separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        character*(*), intent(in):: value
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_string, &
                                     TRIM(ADJUSTL(value))//char(0), &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_string
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_real(attribute, io, name, &
                                                     value, variable_name, &
                                                     separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real, intent(in):: value
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_real, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_real
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_dp(attribute, io, name, &
                                                   value, variable_name, &
                                                   separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(kind=8), intent(in):: value
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_dp, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_dp
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_integer1(attribute, io, name, &
                                                         value, variable_name, &
                                                         separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=1), intent(in):: value
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer1, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer1
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_integer2(attribute, io, name, &
                                                         value, variable_name, &
                                                         separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=2), intent(in):: value
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer2, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer2
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_integer4(attribute, io, name, &
                                                         value, variable_name, &
                                                         separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=4), intent(in):: value
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer4, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer4
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_integer8(attribute, io, name, &
                                                         value, variable_name, &
                                                         separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=8), intent(in):: value
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer8, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer8
            attribute%length = 1
        end if
    end subroutine

    ! 1D
    subroutine adios2_define_variable_attribute_string_1d(attribute, io, name, &
                                                          data, length, &
                                                          variable_name, &
                                                          separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        character*(*), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        ! local data with zero terminated character
        character(len=adios2_string_array_element_max_size), &
            dimension(length):: data_null_terminated

        integer :: i

        do i = 1, length
            data_null_terminated(i) = TRIM(ADJUSTL(data(i)))//char(0)
        end do

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_string, &
                                           data_null_terminated, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_string
            attribute%length = length
        end if

    end subroutine

    subroutine adios2_define_variable_attribute_real_1d(attribute, io, &
                                                        name, data, &
                                                        length, &
                                                        variable_name, &
                                                        separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real, dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_real, data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_real
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_dp_1d(attribute, io, &
                                                      name, data, &
                                                      length, &
                                                      variable_name, &
                                                      separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(kind=8), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_dp, data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_dp
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_integer1_1d(attribute, io, &
                                                            name, data, &
                                                            length, &
                                                            variable_name, &
                                                            separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=1), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer1, &
                                           data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer1
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_integer2_1d(attribute, io, &
                                                            name, data, &
                                                            length, &
                                                            variable_name, &
                                                            separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=2), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer2, &
                                           data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer2
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_integer4_1d(attribute, io, &
                                                            name, data, &
                                                            length, &
                                                            variable_name, &
                                                            separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=4), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer4, data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer4
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_define_variable_attribute_integer8_1d(attribute, io, &
                                                            name, data, &
                                                            length, &
                                                            variable_name, &
                                                            separator, ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=8), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        character*(*), intent(in):: separator
        integer, intent(out) :: ierr

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer8, &
                                           data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer8
            attribute%length = length
        end if
    end subroutine

    ! Variable attribute default separator
    subroutine adios2_def_var_att_sep_string(attribute, io, name, &
                                             value, variable_name, &
                                             ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        character*(*), intent(in):: value
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_string, &
                                     TRIM(ADJUSTL(value))//char(0), &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_string
            attribute%length = 1
        end if

    end subroutine

    subroutine adios2_def_var_att_sep_real(attribute, io, name, &
                                           value, variable_name, &
                                           ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real, intent(in):: value
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_real, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_real
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_dp(attribute, io, name, &
                                         value, variable_name, &
                                         ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(kind=8), intent(in):: value
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_dp, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_dp
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_integer1(attribute, io, name, &
                                               value, variable_name, &
                                               ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=1), intent(in):: value
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer1, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer1
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_integer2(attribute, io, name, &
                                               value, variable_name, &
                                               ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=2), intent(in):: value
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer2, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer2
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_integer4(attribute, io, name, &
                                               value, variable_name, &
                                               ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=4), intent(in):: value
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer4, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer4
            attribute%length = 1
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_integer8(attribute, io, name, &
                                               value, variable_name, &
                                               ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=8), intent(in):: value
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_f2c(attribute%f2c, io%f2c, &
                                     TRIM(ADJUSTL(name))//char(0), &
                                     adios2_type_integer8, value, &
                                     TRIM(ADJUSTL(variable_name))//char(0), &
                                     TRIM(ADJUSTL(separator))//char(0), &
                                     ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .true.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer8
            attribute%length = 1
        end if
    end subroutine

    ! 1D
    subroutine adios2_def_var_att_sep_string_1d(attribute, io, name, &
                                                data, length, &
                                                variable_name, &
                                                ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        character*(*), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        ! local data with zero terminated character
        character(len=adios2_string_array_element_max_size), &
            dimension(length):: data_null_terminated
        integer :: i

        separator = '/'

        do i = 1, length
            data_null_terminated(i) = TRIM(ADJUSTL(data(i)))//char(0)
        end do

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_string, &
                                           data_null_terminated, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_string
            attribute%length = length
        end if

    end subroutine

    subroutine adios2_def_var_att_sep_real_1d(attribute, io, &
                                              name, data, &
                                              length, &
                                              variable_name, &
                                              ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real, dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_real, data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_real
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_dp_1d(attribute, io, &
                                            name, data, &
                                            length, &
                                            variable_name, &
                                            ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        real(kind=8), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_dp, data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_dp
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_integer1_1d(attribute, io, &
                                                  name, data, &
                                                  length, &
                                                  variable_name, &
                                                  ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=1), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer1, &
                                           data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer1
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_integer2_1d(attribute, io, &
                                                  name, data, &
                                                  length, &
                                                  variable_name, &
                                                  ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=2), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer2, &
                                           data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer2
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_integer4_1d(attribute, io, &
                                                  name, data, &
                                                  length, &
                                                  variable_name, &
                                                  ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=4), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer4, data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer4
            attribute%length = length
        end if
    end subroutine

    subroutine adios2_def_var_att_sep_integer8_1d(attribute, io, &
                                                  name, data, &
                                                  length, &
                                                  variable_name, &
                                                  ierr)
        type(adios2_attribute), intent(out) :: attribute
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer(kind=8), dimension(:), intent(in):: data
        integer, intent(in) :: length
        character*(*), intent(in):: variable_name
        integer, intent(out) :: ierr
        !local
        character(len=1):: separator
        separator = '/'

        call adios2_define_vattr_array_f2c(attribute%f2c, io%f2c, &
                                           TRIM(ADJUSTL(name))//char(0), &
                                           adios2_type_integer8, &
                                           data, length, &
                                           TRIM(ADJUSTL(variable_name))//char(0), &
                                           TRIM(ADJUSTL(separator))//char(0), &
                                           ierr)
        if (ierr == 0) then
            attribute%valid = .true.
            attribute%is_value = .false.
            attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(name)
            attribute%type = adios2_type_integer8
            attribute%length = length
        end if
    end subroutine

end module
