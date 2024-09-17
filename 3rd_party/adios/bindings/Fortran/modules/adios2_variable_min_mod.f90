!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_variable_min_mod.f90 : ADIOS2 Fortran bindings for overloaded
!                                adios2_variable_min subroutine
!   Created on: Nov 15, 2018
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_variable_min_mod
    use adios2_parameters_mod
    use adios2_variable_mod
    implicit none

    interface adios2_variable_min

        module procedure adios2_variable_min_real
        module procedure adios2_variable_min_dp
        module procedure adios2_variable_min_complex
        module procedure adios2_variable_min_complex_dp
        module procedure adios2_variable_min_integer1
        module procedure adios2_variable_min_integer2
        module procedure adios2_variable_min_integer4
        module procedure adios2_variable_min_integer8

    end interface
    external adios2_variable_min_f2c

contains

    subroutine adios2_variable_min_real(minimum, variable, ierr)
        real, intent(out) :: minimum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_real, &
                                        'variable_min', ierr)
        if (ierr == 0) then
            call adios2_variable_min_f2c(minimum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_min_dp(minimum, variable, ierr)
        real(kind=8), intent(out) :: minimum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_dp, &
                                        'variable_min', ierr)
        if (ierr == 0) then
            call adios2_variable_min_f2c(minimum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_min_complex(minimum, variable, ierr)
        complex, intent(out) :: minimum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_complex, &
                                        'variable_min', ierr)
        if (ierr == 0) then
            call adios2_variable_min_f2c(minimum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_min_complex_dp(minimum, variable, ierr)
        complex(kind=8), intent(out) :: minimum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_complex_dp, &
                                        'variable_min', ierr)
        if (ierr == 0) then
            call adios2_variable_min_f2c(minimum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_min_integer1(minimum, variable, ierr)
        integer(kind=1), intent(out) :: minimum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_integer1, &
                                        'variable_min', ierr)
        if (ierr == 0) then
            call adios2_variable_min_f2c(minimum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_min_integer2(minimum, variable, ierr)
        integer(kind=2), intent(out) :: minimum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_integer2, &
                                        'variable_min', ierr)
        if (ierr == 0) then
            call adios2_variable_min_f2c(minimum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_min_integer4(minimum, variable, ierr)
        integer(kind=4), intent(out) :: minimum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_integer4, &
                                        'variable_min', ierr)
        if (ierr == 0) then
            call adios2_variable_min_f2c(minimum, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_min_integer8(minimum, variable, ierr)
        integer(kind=8), intent(out) :: minimum
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_check_type(variable, adios2_type_integer8, &
                                        'variable_min', ierr)
        if (ierr == 0) then
            call adios2_variable_min_f2c(minimum, variable%f2c, ierr)
        end if

    end subroutine

end module
