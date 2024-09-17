!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_io_mod.f90 : ADIOS2 Fortran bindings for IO class
!
!   Created on: Mar 13, 2017
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_io_mod
   use adios2_io_open_mod
   use adios2_io_define_variable_mod
   use adios2_io_define_derived_variable_mod
   use adios2_io_define_attribute_mod
   use adios2_variable_mod
   implicit none

   external adios2_add_transport_f2c
   external adios2_attribute_is_value_f2c
   external adios2_attribute_length_f2c
   external adios2_attribute_type_f2c
   external adios2_clear_parameters_f2c
   external adios2_flush_all_engines_f2c
   external adios2_get_parameter_f2c
   external adios2_get_parameter_length_f2c
   external adios2_in_config_file_f2c
   external adios2_available_variables_f2c
   external adios2_available_attributes_f2c
   external adios2_retrieve_namelist_f2c
   external adios2_inquire_attribute_f2c
   external adios2_inquire_variable_attribute_f2c
   external adios2_inquire_variable_f2c
   external adios2_io_engine_type_f2c
   external adios2_io_engine_type_length_f2c
   external adios2_remove_all_attributes_f2c
   external adios2_remove_all_variables_f2c
   external adios2_remove_attribute_f2c
   external adios2_remove_variable_f2c
   external adios2_set_engine_f2c
   external adios2_set_parameter_f2c
   external adios2_set_parameters_f2c
   external adios2_set_transport_parameter_f2c

contains

   subroutine adios2_in_config_file(result, io, ierr)
      logical, intent(out):: result
      type(adios2_io), intent(in):: io
      integer, intent(out):: ierr

      ! local
      integer resultInt

      call adios2_in_config_file_f2c(resultInt, io, ierr)
      if(resultInt == 0) then
         result = .false.
      else
         result = .true.
      end if

   end subroutine

   subroutine adios2_io_engine_type(type, io, ierr)
      character(len=:), allocatable, intent(out) :: type
      type(adios2_io), intent(in) :: io
      integer, intent(out) :: ierr

      !local
      integer :: length

      if (allocated(type)) deallocate (type)

      call adios2_io_engine_type_length_f2c(length, io%f2c, ierr)
      if (ierr == 0) then
         allocate (character(length) :: type)
         call adios2_io_engine_type_f2c(type, io%f2c, ierr)
      end if

   end subroutine

   subroutine adios2_set_engine(io, engine_type, ierr)
      type(adios2_io), intent(inout) :: io
      character*(*), intent(in) :: engine_type
      integer, intent(out) :: ierr

      call adios2_set_engine_f2c(io%f2c, &
         TRIM(ADJUSTL(engine_type))//char(0), ierr)

      if( ierr == 0 ) io%engine_type = engine_type

   end subroutine

   subroutine adios2_set_parameters(io, parameters, ierr)
      type(adios2_io), intent(in) :: io
      character*(*), intent(in) :: parameters
      integer, intent(out) :: ierr

      call adios2_set_parameters_f2c(io%f2c, &
         TRIM(ADJUSTL(parameters))//char(0), &
         ierr)
   end subroutine

   subroutine adios2_set_parameter(io, key, value, ierr)
      type(adios2_io), intent(in) :: io
      character*(*), intent(in) :: key
      character*(*), intent(in) :: value
      integer, intent(out) :: ierr

      call adios2_set_parameter_f2c(io%f2c, TRIM(ADJUSTL(key))//char(0), &
         TRIM(ADJUSTL(value))//char(0), ierr)
   end subroutine

   subroutine adios2_get_parameter(value, io, key, ierr)
      character(len=:), allocatable, intent(out) :: value
      type(adios2_io), intent(in) :: io
      character*(*), intent(in) :: key
      integer, intent(out) :: ierr

      !local
      integer :: length

      if (allocated(value)) deallocate (value)

      call adios2_get_parameter_length_f2c(length, io%f2c, &
         TRIM(ADJUSTL(key))//char(0), ierr)
      if (ierr == 0) then
         allocate (character(length) :: value)
         call adios2_get_parameter_f2c(value, io%f2c, &
            TRIM(ADJUSTL(key))//char(0), ierr)
      end if
   end subroutine

   subroutine adios2_clear_parameters(io, ierr)
      type(adios2_io), intent(in) :: io
      integer, intent(out) :: ierr
      call adios2_clear_parameters_f2c(io%f2c, ierr)
   end subroutine

   subroutine adios2_add_transport(transport_index, io, type, ierr)
      integer, intent(out):: transport_index
      type(adios2_io), intent(in) :: io
      character*(*), intent(in) :: type
      integer, intent(out) :: ierr

      call adios2_add_transport_f2c(transport_index, io%f2c, &
         TRIM(ADJUSTL(type))//char(0), ierr)

   end subroutine

   subroutine adios2_set_transport_parameter(io, transport_index, key, value, &
      ierr)
      type(adios2_io), intent(in):: io
      integer, intent(in):: transport_index
      character*(*), intent(in) :: key
      character*(*), intent(in) :: value
      integer, intent(out):: ierr

      call adios2_set_transport_parameter_f2c(io%f2c, transport_index, &
         TRIM(ADJUSTL(key))//char(0), &
         TRIM(ADJUSTL(value))//char(0), &
         ierr)
   end subroutine


   subroutine adios2_available_variables(io, namestruct, ierr)
      type(adios2_io), intent(in) :: io
      type(adios2_namestruct), intent(out) :: namestruct
      integer, intent(out) :: ierr

      call adios2_available_variables_f2c(io%f2c, namestruct%f2c, &
         namestruct%count,  namestruct%max_name_len, ierr)
      if (ierr == 0) then
         namestruct%valid = .true.
      endif
   end subroutine

   subroutine adios2_retrieve_names(namestruct, namelist, ierr)
      type(adios2_namestruct), intent(inout) :: namestruct
      character(*), dimension(*), intent(inout) :: namelist
      integer, intent(out) :: ierr

      if (namestruct%valid .and. namestruct%f2c > 0_8) then
         call adios2_retrieve_namelist_f2c(namestruct%f2c, namelist, ierr)
      else
         write(*,*) "ADIOS2 Fortran ERROR: invalid namestruct when calling adios2_retrieve_names()"
      endif
      namestruct%valid = .false.
   end subroutine

   !
   ! F2008 implementation that allows for allocating a character array inside
   !
   ! subroutine adios2_available_variables(io, nvars, varnamelist, ierr)
   !     type(adios2_io), intent(in) :: io
   !     integer, intent(out) :: nvars
   !     character(len=:), dimension(:), allocatable, intent(out) :: varnamelist
   !     integer, intent(out) :: ierr

   !     integer(kind=8):: namestruct
   !     integer :: count, max_name_len

   !     call adios2_available_variables_f2c(io%f2c, namestruct, count, &
   !                     max_name_len, ierr)
   !     if (ierr == 0) then
   !         allocate(character(len=max_name_len) :: varnamelist(count))
   !     endif

   !     call adios2_retrieve_variable_names_f2c(namestruct, varnamelist, ierr)
   !     nvars = count
   ! end subroutine


   subroutine adios2_inquire_variable(variable, io, name, ierr)
      type(adios2_variable), intent(out) :: variable
      type(adios2_io), intent(in) :: io
      character*(*), intent(in) :: name
      integer, intent(out) :: ierr

      call adios2_inquire_variable_f2c(variable%f2c, io%f2c, &
         TRIM(ADJUSTL(name))//char(0), ierr)

      if(variable%f2c > 0_8) then
         variable%valid = .true.
         variable%name = name
         call adios2_variable_type(variable%type, variable, ierr)
         call adios2_variable_ndims(variable%ndims, variable, ierr)
      else
         variable%valid = .false.
         variable%name = ''
         variable%type = adios2_type_unknown
         variable%ndims = -1
      end if

   end subroutine

   subroutine adios2_remove_variable(result, io, name, ierr)
      type(adios2_io), intent(in) :: io
      character*(*), intent(in) :: name
      logical, intent(out) :: result
      integer, intent(out) :: ierr
      ! Local
      type(adios2_variable):: variable
      integer:: resultInt

      call adios2_inquire_variable(variable, io, name, ierr)
      if( variable%valid ) then
         call adios2_remove_variable_f2c(resultInt, io%f2c, &
            TRIM(ADJUSTL(name))//char(0), ierr)
         if( resultInt == 1) then
            result = .true.
         else
            result = .false.
         end if
      end if

   end subroutine


   subroutine adios2_remove_all_variables(io, ierr)
      type(adios2_io), intent(in) :: io
      integer, intent(out) :: ierr

      call adios2_remove_all_variables_f2c(io%f2c, ierr)

   end subroutine

   subroutine adios2_available_attributes(io, namestruct, ierr)
      type(adios2_io), intent(in) :: io
      type(adios2_namestruct), intent(out) :: namestruct
      integer, intent(out) :: ierr

      call adios2_available_attributes_f2c(io%f2c, namestruct%f2c, &
         namestruct%count,  namestruct%max_name_len, ierr)
      if (ierr == 0) then
         namestruct%valid = .true.
      endif
   end subroutine

   ! subroutine adios2_available_attributes(io, nattrs, attrnamelist, ierr)
   !     type(adios2_io), intent(in) :: io
   !     integer, intent(out) :: nattrs
   !     character(len=:), dimension(:), allocatable, intent(out) :: attrnamelist
   !     integer, intent(out) :: ierr

   !     integer(kind=8):: namestruct
   !     integer :: count, max_name_len

   !     call adios2_available_attributes_f2c(io%f2c, namestruct, count, &
   !                     max_name_len, ierr)
   !     if (ierr == 0) then
   !         allocate(character(len=max_name_len) :: attrnamelist(count))
   !     endif

   !     call adios2_retrieve_attribute_names_f2c(namestruct, count, &
   !                     max_name_len, attrnamelist, ierr)
   !     nattrs = count
   ! end subroutine

   subroutine adios2_inquire_attribute(attribute, io, name, ierr)
      type(adios2_attribute), intent(out) :: attribute
      type(adios2_io), intent(in) :: io
      character*(*), intent(in) :: name
      integer, intent(out) :: ierr
      !local
      integer:: is_valueInt

      call adios2_inquire_attribute_f2c(attribute%f2c, io%f2c, &
         TRIM(ADJUSTL(name))//char(0), ierr)

      if(attribute%f2c > 0_8) then
         attribute%valid = .true.
         attribute%name = name
         call adios2_attribute_type_f2c(attribute%type, attribute%f2c, ierr)
         call adios2_attribute_length_f2c(attribute%length, attribute%f2c, &
            ierr)
         call adios2_attribute_is_value_f2c(is_valueInt, attribute%f2c, ierr)

         if(is_valueInt == 0) then
            attribute%is_value = .false.
         else
            attribute%is_value = .true.
         end if

      else
         attribute%valid = .false.
         attribute%name = ''
         attribute%type = adios2_type_unknown
         attribute%length = 0
      end if

   end subroutine

   subroutine adios2_inquire_variable_attribute(attribute, io, attribute_name, variable_name, separator, ierr)
      type(adios2_attribute), intent(out) :: attribute
      type(adios2_io), intent(in)         :: io
      character*(*), intent(in)           :: attribute_name
      character*(*), intent(in)           :: variable_name
      character*(*), intent(in)           :: separator
      integer, intent(out)                :: ierr
      !local
      integer:: is_valueInt

      call adios2_inquire_variable_attribute_f2c(attribute%f2c, io%f2c, &
         TRIM(ADJUSTL(attribute_name))//char(0), &
         TRIM(ADJUSTL(variable_name))//char(0), &
         TRIM(ADJUSTL(separator))//char(0), &
         ierr)

      if(attribute%f2c > 0_8) then
         attribute%valid = .true.
         attribute%name = TRIM(variable_name)//TRIM(separator)//TRIM(attribute_name)
         call adios2_attribute_type_f2c(attribute%type, attribute%f2c, ierr)
         call adios2_attribute_length_f2c(attribute%length, attribute%f2c, &
            ierr)
         call adios2_attribute_is_value_f2c(is_valueInt, attribute%f2c, ierr)

         if(is_valueInt == 0) then
            attribute%is_value = .false.
         else
            attribute%is_value = .true.
         end if

      else
         attribute%valid = .false.
         attribute%name = ''
         attribute%type = adios2_type_unknown
         attribute%length = 0
      end if

   end subroutine


   subroutine adios2_remove_attribute(result, io, name, ierr)
      type(adios2_io), intent(in) :: io
      character*(*), intent(in) :: name
      logical, intent(out) :: result
      integer, intent(out) :: ierr

      ! Local
      integer :: resultInt

      call adios2_remove_attribute_f2c(resultInt, io%f2c, &
         TRIM(ADJUSTL(name))//char(0), ierr)
      if( resultInt == 1) then
         result = .true.
      else
         result = .false.
      end if

   end subroutine


   subroutine adios2_remove_all_attributes(io, ierr)
      type(adios2_io), intent(in) :: io
      integer, intent(out) :: ierr

      call adios2_remove_all_attributes_f2c(io%f2c, ierr)

   end subroutine


   subroutine adios2_flush_all_engines(io, ierr)
      type(adios2_io), intent(in) :: io
      integer, intent(out) :: ierr

      call adios2_flush_all_engines_f2c(io%f2c, ierr)

   end subroutine

end module
