!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_attribute_mod.f90 : ADIOS2 Fortran bindings for
!                             type(adios2_attribute) handler subroutines
!   Created on: Dec 10, 2018
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_attribute_mod
    use adios2_parameters_mod
    implicit none

    contains

    subroutine adios2_attribute_name(name, attribute, ierr)
        character(len=:), allocatable, intent(out) :: name
        type(adios2_attribute), intent(in) :: attribute
        integer, intent(out) :: ierr
        external adios2_attribute_name_length_f2c
        external adios2_attribute_name_f2c

        !local
        integer :: length

        if (allocated(name)) deallocate (name)

        call adios2_attribute_name_length_f2c(length, attribute%f2c, ierr)
        if (ierr == 0) then
            allocate (character(length) :: name)
            call adios2_attribute_name_f2c(name, attribute%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_attribute_check_type(attribute, adios2_type, hint, ierr)
        type(adios2_attribute), intent(in):: attribute
        integer, intent(in):: adios2_type
        character*(*), intent(in):: hint
        integer, intent(out):: ierr

        if (attribute%type /= adios2_type) then
            write (0, *) 'ERROR: adios2 attribute ', TRIM(attribute%name)//char(0), &
                ' type mismatch, in call to adios2_', TRIM(hint)//char(0)

            ierr = adios2_error_invalid_argument
        end if

    end subroutine

end module
