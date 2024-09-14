!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_variable_mod.f90 : ADIOS2 Fortran bindings for Variable class
!
!   Created on: Mar 13, 2017
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_variable_mod
    use adios2_functions_mod
    implicit none
    external adios2_add_operation_f2c
    external adios2_remove_operations_f2c
    external adios2_set_block_selection_f2c
    external adios2_set_memory_selection_f2c
    external adios2_set_operation_parameter_f2c
    external adios2_set_selection_f2c
    external adios2_set_shape_f2c
    external adios2_set_step_selection_f2c
    external adios2_variable_name_f2c
    external adios2_variable_name_length_f2c
    external adios2_variable_ndims_f2c
    external adios2_variable_shape_f2c
    external adios2_variable_steps_f2c
    external adios2_variable_type_f2c
contains

    subroutine adios2_variable_name(name, variable, ierr)
        character(len=:), allocatable, intent(out) :: name
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        !local
        integer :: length

        if (allocated(name)) deallocate (name)

        call adios2_variable_name_length_f2c(length, variable%f2c, ierr)
        if (ierr == 0) then
            allocate (character(length) :: name)
            call adios2_variable_name_f2c(name, variable%f2c, ierr)
        end if

    end subroutine

    subroutine adios2_variable_type(type, variable, ierr)
        integer, intent(out) :: type
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr
        ! local
        integer :: c_type

        call adios2_variable_type_f2c(c_type, variable%f2c, ierr)
        call adios2_TypeC2F(c_type, type)

    end subroutine

    subroutine adios2_variable_ndims(ndims, variable, ierr)
        integer, intent(out) :: ndims
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_ndims_f2c(ndims, variable%f2c, ierr)
    end subroutine

    subroutine adios2_variable_shape(shape_dims, ndims, variable, ierr)
        integer(kind=8), dimension(:), allocatable, intent(out) :: shape_dims
        integer, intent(out) :: ndims
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_ndims_f2c(ndims, variable%f2c, ierr)
        allocate (shape_dims(ndims))
        call adios2_variable_shape_f2c(shape_dims, variable%f2c, ierr)

    end subroutine

    subroutine adios2_variable_steps(steps, variable, ierr)
        integer(kind=8), intent(out) :: steps
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_variable_steps_f2c(steps, variable%f2c, ierr)

    end subroutine

    subroutine adios2_set_memory_space(variable, mem, ierr)
        type(adios2_variable), intent(in) :: variable
        integer, intent(in) :: mem
        integer, intent(out) :: ierr

        call adios2_set_memory_space_f2c(variable%f2c, mem, ierr)
    end subroutine

    subroutine adios2_get_memory_space(mem, variable, ierr)
        integer, intent(out) :: mem
        type(adios2_variable), intent(in) :: variable
        integer, intent(out) :: ierr

        call adios2_get_memory_space_f2c(mem, variable%f2c, ierr)
    end subroutine

    subroutine adios2_set_shape(variable, ndims, shape_dims, ierr)
        type(adios2_variable), intent(in) :: variable
        integer, intent(in) :: ndims
        integer(kind=8), dimension(:), intent(in) :: shape_dims
        integer, intent(out) :: ierr

        call adios2_set_shape_f2c(variable%f2c, ndims, shape_dims, ierr)
    end subroutine

    subroutine adios2_set_block_selection(variable, block_id, ierr)
        type(adios2_variable), intent(in) :: variable
        integer(kind=8), intent(in) :: block_id
        integer, intent(out) :: ierr

        call adios2_set_block_selection_f2c(variable%f2c, block_id, ierr)
    end subroutine

    subroutine adios2_set_selection(variable, ndims, start_dims, count_dims, &
                                    ierr)
        type(adios2_variable), intent(in) :: variable
        integer, intent(in) :: ndims
        integer(kind=8), dimension(:), intent(in) :: start_dims
        integer(kind=8), dimension(:), intent(in) :: count_dims
        integer, intent(out) :: ierr

        call adios2_set_selection_f2c(variable%f2c, ndims, start_dims, &
                                      count_dims, ierr)
    end subroutine

    subroutine adios2_set_memory_selection(variable, ndims, &
                                           memory_start_dims, &
                                           memory_count_dims, &
                                           ierr)
        type(adios2_variable), intent(in) :: variable
        integer, intent(in) :: ndims
        integer(kind=8), dimension(:), intent(in) :: memory_start_dims
        integer(kind=8), dimension(:), intent(in) :: memory_count_dims
        integer, intent(out) :: ierr

        call adios2_set_memory_selection_f2c(variable%f2c, ndims, &
                                             memory_start_dims, &
                                             memory_count_dims, ierr)
    end subroutine

    subroutine adios2_set_step_selection(variable, step_start, step_count, ierr)
        type(adios2_variable), intent(in) :: variable
        integer(kind=8), intent(in) :: step_start
        integer(kind=8), intent(in) :: step_count
        integer, intent(out) :: ierr

        call adios2_set_step_selection_f2c(variable%f2c, step_start, &
                                           step_count, ierr)
    end subroutine


    subroutine adios2_variable_check_type(variable, adios2_type, hint, ierr)
        type(adios2_variable), intent(in):: variable
        integer, intent(in):: adios2_type
        character*(*), intent(in):: hint
        integer, intent(out):: ierr

        if( variable%type /= adios2_type ) then
            write(0,*) 'ERROR: adios2 variable ', TRIM(variable%name)//char(0), &
                       ' type mismatch, in call to adios2_', TRIM(hint)//char(0), &
                       'variable type: ', variable%type, ' expected type: ', &
                       adios2_type

            ierr = adios2_error_invalid_argument
        else
            ierr = adios2_error_none
        end if

    end subroutine

    subroutine adios2_add_operation(operation_index, variable, op, key, value, &
                                    ierr)
        integer, intent(out):: operation_index
        type(adios2_variable), intent(in):: variable
        type(adios2_operator), intent(in):: op
        character*(*), intent(in):: key
        character*(*), intent(in):: value
        integer, intent(out):: ierr

        call adios2_add_operation_f2c(operation_index, variable%f2c, op%f2c, &
                                      TRIM(ADJUSTL(key))//char(0), &
                                      TRIM(ADJUSTL(value))//char(0), ierr)
    end subroutine


    subroutine adios2_set_operation_parameter(variable, operation_index, key, &
                                              value, ierr)
        type(adios2_variable), intent(in):: variable
        integer, intent(in):: operation_index
        character*(*), intent(in):: key
        character*(*), intent(in):: value
        integer, intent(out):: ierr

        call adios2_set_operation_parameter_f2c(variable%f2c, operation_index, &
                                                TRIM(ADJUSTL(key))//char(0), &
                                                TRIM(ADJUSTL(value))//char(0), &
                                                ierr)
    end subroutine

    subroutine adios2_remove_operations(variable, ierr)
        type(adios2_variable), intent(in):: variable
        integer, intent(out):: ierr

        call adios2_remove_operations_f2c(variable%f2c, ierr)
    end subroutine

end module
