!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_io_define_variable_mod.f90 : ADIOS2 Fortran bindings for IO class
!  overloaded (C++ template) function adios2_define_variable
!
!   Created on: Mar 13, 2017
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_io_define_variable_mod
    use adios2_parameters_mod
    use adios2_functions_mod
    implicit none

    interface adios2_define_variable

        module procedure adios2_define_variable_value
        module procedure adios2_define_variable_array

    end interface

    external adios2_define_global_variable_f2c
    external adios2_define_variable_f2c

contains

    subroutine adios2_define_variable_value(variable, io, name, adios2_type, &
                                            ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer, intent(in):: adios2_type
        integer, intent(out) :: ierr

        call adios2_define_global_variable_f2c(variable%f2c, io%f2c, &
                                               TRIM(ADJUSTL(name))//char(0), &
                                               adios2_type, ierr)
        if( ierr == 0 ) then
            variable%valid = .true.
            variable%name = name
            variable%type = adios2_type
            variable%ndims = 1
        end if

    end subroutine

    subroutine adios2_define_variable_array(variable, io, name, adios2_type, &
                                            ndims, shape_dims, start_dims, &
                                            count_dims, is_constant_dims, ierr)
        type(adios2_variable), intent(out) :: variable
        type(adios2_io), intent(in) :: io
        character*(*), intent(in) :: name
        integer, intent(in) :: adios2_type
        integer, intent(in) :: ndims
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        logical, intent(in) :: is_constant_dims
        integer, intent(out) :: ierr

        integer is_constant_dims_int
        is_constant_dims_int = adios2_LogicalToInt(is_constant_dims)

        call adios2_define_variable_f2c(variable%f2c, io%f2c, &
                                        TRIM(ADJUSTL(name))//char(0), &
                                        adios2_type, ndims, &
                                        shape_dims, start_dims, count_dims, &
                                        is_constant_dims_int, ierr)

        if( ierr == 0 ) then
            variable%valid = .true.
            variable%name = name
            variable%type = adios2_type
            variable%ndims = ndims
        end if

    end subroutine

end module
