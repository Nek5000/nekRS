!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_adios_mod.f90 : ADIOS2 Fortran bindings for the ADIOS component
!
!   Created on: Aug 22, 2017
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_adios_mod
    use adios2_adios_init_mod
    implicit none

    external adios2_declare_io_f2c
    external adios2_io_engine_type_length_f2c
    external adios2_io_engine_type_f2c
    external adios2_at_io_f2c
    external adios2_define_operator_f2c
    external adios2_inquire_operator_f2c
    external adios2_operator_type_length_f2c
    external adios2_operator_type_f2c
    external adios2_flush_all_f2c
    external adios2_remove_io_f2c
    external adios2_remove_all_ios_f2c
    external adios2_finalize_f2c
    external adios2_enter_computation_block_f2c
    external adios2_exit_computation_block_f2c
contains

    subroutine adios2_declare_io(io, adios, io_name, ierr)
        type(adios2_io), intent(inout) :: io
        type(adios2_adios), intent(in) :: adios
        character*(*), intent(in)  :: io_name
        integer, intent(out) :: ierr
        !local
        integer:: length

        call adios2_declare_io_f2c(io%f2c, adios%f2c, &
                                   TRIM(ADJUSTL(io_name))//char(0), ierr)
        if( ierr == 0 ) then
            io%valid = .true.
            call adios2_io_engine_type_length_f2c(length, io%f2c, ierr)
            if (length > 15) stop 'adios2_declare_io: engine_type too long!'
            call adios2_io_engine_type_f2c(io%engine_type, io%f2c, ierr)
            io%engine_type = io%engine_type(1:length)
        end if

    end subroutine

    subroutine adios2_at_io(io, adios, io_name, ierr)
        type(adios2_io), intent(out) :: io
        type(adios2_adios), intent(in) :: adios
        character*(*), intent(in)  :: io_name
        integer, intent(out) :: ierr
        !local
        integer:: length

        call adios2_at_io_f2c(io%f2c, adios%f2c, &
                              TRIM(ADJUSTL(io_name))//char(0), ierr)
        if( ierr == 0 ) then
            io%valid = .true.
            call adios2_io_engine_type_length_f2c(length, io%f2c, ierr)
            if (length > 15) stop 'adios2_at_io: engine_type too long!'
            call adios2_io_engine_type_f2c(io%engine_type, io%f2c, ierr)
            io%engine_type = io%engine_type(1:length)
        end if

    end subroutine

    subroutine adios2_define_operator(op, adios, op_name, op_type, ierr)
        type(adios2_operator), intent(out) :: op
        type(adios2_adios), intent(in) :: adios
        character*(*), intent(in)  :: op_name
        character*(*), intent(in)  :: op_type
        integer, intent(out) :: ierr

        call adios2_define_operator_f2c(op%f2c, adios%f2c, &
                                        TRIM(ADJUSTL(op_name))//char(0), &
                                        TRIM(ADJUSTL(op_type))//char(0), ierr)
        if(ierr == 0) then
            op%valid = .true.
            op%name = op_name
            op%type = op_type
        end if

    end subroutine

    subroutine adios2_inquire_operator(op, adios, op_name, ierr)
        type(adios2_operator), intent(out) :: op
        type(adios2_adios), intent(in) :: adios
        character*(*), intent(in)  :: op_name
        integer, intent(out) :: ierr
        !local
        integer :: length

        call adios2_inquire_operator_f2c(op%f2c, adios%f2c, &
                                         TRIM(ADJUSTL(op_name))//char(0), ierr)

        if(ierr == adios2_found) then
            op%valid = .true.
            op%name = op_name

            call adios2_operator_type_length_f2c(length, op%f2c, ierr)
            if (length > 64) stop 'adios2_inquire_operator: operator_type too long!'
            call adios2_operator_type_f2c(op%type, op%f2c, ierr)
            op%type = op%type(1:length)
        else
            op%valid = .false.
            op%name = ''
            op%type = ''
        end if

    end subroutine

    subroutine adios2_flush_all(adios, ierr)
        type(adios2_adios), intent(in) :: adios
        integer, intent(out) :: ierr

        call adios2_flush_all_f2c(adios%f2c, ierr)

    end subroutine

    subroutine adios2_remove_io(result, adios, name, ierr)
        logical, intent(out):: result
        type(adios2_adios), intent(in) :: adios
        character*(*), intent(in):: name
        integer, intent(out) :: ierr
        ! local
        integer resultInt

        call adios2_remove_io_f2c(resultInt, adios, &
                                  TRIM(ADJUSTL(name))//char(0), ierr)
        if(resultInt == 0) then
            result = .false.
        else
            result = .true.
        end if

    end subroutine

    subroutine adios2_remove_all_ios(adios, ierr)
        type(adios2_adios), intent(in) :: adios
        integer, intent(out) :: ierr

        call adios2_remove_all_ios_f2c(adios, ierr)

    end subroutine

    subroutine adios2_finalize(adios, ierr)
        type(adios2_adios), intent(inout) :: adios
        integer, intent(out) :: ierr

        call adios2_finalize_f2c(adios%f2c, ierr)
        adios%valid = .false.

    end subroutine


    subroutine adios2_enter_computation_block(adios, ierr)
        type(adios2_adios), intent(in) :: adios
        integer, intent(out) :: ierr
        call adios2_enter_computation_block_f2c(adios, ierr)
    end subroutine

    subroutine adios2_exit_computation_block(adios, ierr)
        type(adios2_adios), intent(in) :: adios
        integer, intent(out) :: ierr
        call adios2_exit_computation_block_f2c(adios, ierr)
    end subroutine


end module
