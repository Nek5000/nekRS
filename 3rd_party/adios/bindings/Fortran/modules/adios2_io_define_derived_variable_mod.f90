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

module adios2_io_define_derived_variable_mod
   use adios2_parameters_mod
   implicit none

   external adios2_define_derived_variable_f2c

contains

   subroutine adios2_define_derived_variable(variable, io, name, expression, adios2_derived_var_type, &
      ierr)
      type(adios2_derived_variable), intent(out) :: variable
      type(adios2_io), intent(in) :: io
      character*(*), intent(in) :: name
      character*(*), intent(in) :: expression
      integer, intent(in):: adios2_derived_var_type
      integer, intent(out) :: ierr

      call adios2_define_derived_variable_f2c(variable%f2c, io%f2c, &
         TRIM(ADJUSTL(name))//char(0), TRIM(ADJUSTL(expression))//char(0), &
         adios2_derived_var_type, ierr)
      if( ierr == 0 ) then
         variable%valid = .true.
         variable%name = name
         variable%type = adios2_derived_var_type
      end if

   end subroutine

end module
