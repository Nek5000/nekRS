!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_mod.f90 : ADIOS2 Fortran bindings central module
!
!   Created on: Mar 13, 2017
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_parameters_mod
   implicit none

   ! Types
   integer, parameter :: adios2_type_unknown = -1

   integer, parameter :: adios2_type_character = 0
   integer, parameter :: adios2_type_string = 0

   integer, parameter :: adios2_type_real4 = 1
   integer, parameter :: adios2_type_real  = 1

   integer, parameter :: adios2_type_real8 = 2
   integer, parameter :: adios2_type_dp    = 2
   integer, parameter :: adios2_type_double_precision = 2

   integer, parameter :: adios2_type_complex4 = 3
   integer, parameter :: adios2_type_complex = 3

   integer, parameter :: adios2_type_complex8 = 4
   integer, parameter :: adios2_type_complex_dp = 4

   integer, parameter :: adios2_type_integer1 = 5
   integer, parameter :: adios2_type_integer2 = 6
   integer, parameter :: adios2_type_integer4 = 7
   integer, parameter :: adios2_type_integer8 = 8

   ! is_constant_dims
   logical, parameter :: adios2_constant_dims = .true.
   logical, parameter :: adios2_variable_dims = .false.

   ! Variable Found or not found, ierr value
   integer, parameter :: adios2_not_found = -1
   integer, parameter :: adios2_found = 0

   ! error
   integer, parameter :: adios2_error_none = 0
   integer, parameter :: adios2_error_invalid_argument = 1
   integer, parameter :: adios2_error_system_error = 2
   integer, parameter :: adios2_error_runtime_error = 3
   integer, parameter :: adios2_error_exception = 4

   ! Mode
   integer, parameter :: adios2_mode_undefined = 0
   integer, parameter :: adios2_mode_write = 1
   integer, parameter :: adios2_mode_read = 2
   integer, parameter :: adios2_mode_append = 3
   integer, parameter :: adios2_mode_readRandomAccess = 6

   integer, parameter :: adios2_mode_deferred = 4
   integer, parameter :: adios2_mode_sync = 5

   integer, parameter :: adios2_memory_space_detect = 0
   integer, parameter :: adios2_memory_space_host = 1
   integer, parameter :: adios2_memory_space_gpu = 2

   ! Step Mode
   integer, parameter :: adios2_step_mode_append = 0
   integer, parameter :: adios2_step_mode_update = 1
   integer, parameter :: adios2_step_mode_read = 2

   ! Step Status
   integer, parameter :: adios2_step_status_other_error = -1
   integer, parameter :: adios2_step_status_ok = 0
   integer, parameter :: adios2_step_status_not_ready = 1
   integer, parameter :: adios2_step_status_end_of_stream = 2

   ! Derived variable type
   integer, parameter :: adios2_derived_var_type_metadata_only = 0
   integer, parameter :: adios2_derived_var_type_expression_string = 1
   integer, parameter :: adios2_derived_var_type_store_data = 2

   !> Fixed size for string array, used in variables and attributes,
   !! must be less or equal than C equivalent in adios2_c_types.h
   integer, parameter :: adios2_string_array_element_max_size = 4096

   integer(kind=8), parameter, dimension(1) :: adios2_null_dims = (/-1/)
   integer(kind=8), parameter :: adios2_local_value_dim = -2

   logical, parameter :: adios2_advance_yes = .true.
   logical, parameter :: adios2_advance_no  = .false.

   ! Low level API handlers
   type adios2_adios
      integer(kind=8):: f2c = 0_8
      logical :: valid = .false.
   end type

   type adios2_io
      integer(kind=8):: f2c = 0_8
      logical :: valid = .false.
      character(len=15):: engine_type = 'BPFile'
   end type

   type adios2_variable
      integer(kind=8):: f2c = 0_8
      logical :: valid = .false.
      character(len=4096):: name = ''
      integer :: type = -1
      integer :: ndims = -1
   end type

   type adios2_derived_variable
      integer(kind=8):: f2c = 0_8
      logical :: valid = .false.
      character(len=4096):: name = ''
      integer :: type = -1
   end type

   type adios2_attribute
      integer(kind=8):: f2c = 0_8
      logical :: valid = .false.
      logical :: is_value = .false.
      character(len=4096):: name = ''
      integer :: type = -1
      integer :: length = -1
   end type

   type adios2_engine
      integer(kind=8):: f2c = 0_8
      logical :: valid = .false.
      character(len=64):: name = ''
      character(len=15):: type = ''
      integer :: mode = adios2_mode_undefined
   end type

   type adios2_operator
      integer(kind=8):: f2c = 0_8
      logical :: valid = .false.
      character(len=64):: name = ''
      character(len=64):: type = ''
   end type

   type adios2_namestruct
      integer(kind=8):: f2c = 0_8
      logical :: valid = .false.
      integer :: count
      integer :: max_name_len
   end type

end module
