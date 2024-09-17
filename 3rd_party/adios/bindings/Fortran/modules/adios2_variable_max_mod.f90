!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_variable_max_mod.f90 : ADIOS2 Fortran bindings for overloaded
!                                adios2_variable_max subroutine
!   Created on: Nov 15, 2018
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_variable_max_mod
    use adios2_parameters_mod
    use adios2_variable_mod
    implicit none

    interface adios2_variable_max

        module procedure adios2_variable_max_real
        module procedure adios2_variable_max_dp
        module procedure adios2_variable_max_complex
        module procedure adios2_variable_max_complex_dp
        module procedure adios2_variable_max_integer1
        module procedure adios2_variable_max_integer2
        module procedure adios2_variable_max_integer4
        module procedure adios2_variable_max_integer8

    end interface

    external adios2_variable_max_f2c
contains

    subroutine adios2_variable_max_real(maximum, variable, ierr)
        real, intent(out) :: maximum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_real, &
                                        'variable_max', ierr)
        if (ierr == 0) then
            call adios2_variable_max_f2c(maximum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_max_dp(maximum, variable, ierr)
        real(kind=8), intent(out) :: maximum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_dp, &
                                        'variable_max', ierr)
        if (ierr == 0) then
            call adios2_variable_max_f2c(maximum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_max_complex(maximum, variable, ierr)
        complex, intent(out) :: maximum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_complex, &
                                        'variable_max', ierr)
        if (ierr == 0) then
            call adios2_variable_max_f2c(maximum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_max_complex_dp(maximum, variable, ierr)
        complex(kind=8), intent(out) :: maximum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_complex_dp, &
                                        'variable_max', ierr)
        if (ierr == 0) then
            call adios2_variable_max_f2c(maximum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_max_integer1(maximum, variable, ierr)
        integer(kind=1), intent(out) :: maximum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_integer1, &
                                        'variable_max', ierr)
        if (ierr == 0) then
            call adios2_variable_max_f2c(maximum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_max_integer2(maximum, variable, ierr)
        integer(kind=2), intent(out) :: maximum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_integer2, &
                                        'variable_max', ierr)
        if (ierr == 0) then
            call adios2_variable_max_f2c(maximum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_max_integer4(maximum, variable, ierr)
        integer(kind=4), intent(out) :: maximum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_integer4, &
                                        'variable_max', ierr)
        if (ierr == 0) then
            call adios2_variable_max_f2c(maximum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_max_integer8(maximum, variable, ierr)
        integer(kind=8), intent(out) :: maximum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_integer8, &
                                        'variable_max', ierr)
        if (ierr == 0) then
            call adios2_variable_max_f2c(maximum, variable%f2c, ierr)
        end if

    end subroutine

end module
