!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_engine_get_mod.f90 : ADIOS2 Fortran bindings for Engine generic
!                              get functions
!
!   Created on: Aug 22, 2017
!       Author: William F Godoy godoywf@ornl.gov
!

module adios2_engine_get_mod
    use adios2_variable_mod
    use adios2_parameters_mod
    implicit none

    interface adios2_get

        ! Single Value
        module procedure adios2_get_string
        module procedure adios2_get_real
        module procedure adios2_get_dp
        module procedure adios2_get_complex
        module procedure adios2_get_complex_dp
        module procedure adios2_get_integer1
        module procedure adios2_get_integer2
        module procedure adios2_get_integer4
        module procedure adios2_get_integer8

        ! 1D Array
        module procedure adios2_get_real_1d
        module procedure adios2_get_dp_1d
        module procedure adios2_get_complex_1d
        module procedure adios2_get_complex_dp_1d
        module procedure adios2_get_integer1_1d
        module procedure adios2_get_integer2_1d
        module procedure adios2_get_integer4_1d
        module procedure adios2_get_integer8_1d

        ! 2D Array
        module procedure adios2_get_real_2d
        module procedure adios2_get_dp_2d
        module procedure adios2_get_complex_2d
        module procedure adios2_get_complex_dp_2d
        module procedure adios2_get_integer1_2d
        module procedure adios2_get_integer2_2d
        module procedure adios2_get_integer4_2d
        module procedure adios2_get_integer8_2d

        ! 3D Array
        module procedure adios2_get_real_3d
        module procedure adios2_get_dp_3d
        module procedure adios2_get_complex_3d
        module procedure adios2_get_complex_dp_3d
        module procedure adios2_get_integer1_3d
        module procedure adios2_get_integer2_3d
        module procedure adios2_get_integer4_3d
        module procedure adios2_get_integer8_3d

        ! 4D Array
        module procedure adios2_get_real_4d
        module procedure adios2_get_dp_4d
        module procedure adios2_get_complex_4d
        module procedure adios2_get_complex_dp_4d
        module procedure adios2_get_integer1_4d
        module procedure adios2_get_integer2_4d
        module procedure adios2_get_integer4_4d
        module procedure adios2_get_integer8_4d

        ! 5D Array
        module procedure adios2_get_real_5d
        module procedure adios2_get_dp_5d
        module procedure adios2_get_complex_5d
        module procedure adios2_get_complex_dp_5d
        module procedure adios2_get_integer1_5d
        module procedure adios2_get_integer2_5d
        module procedure adios2_get_integer4_5d
        module procedure adios2_get_integer8_5d

        ! 6D Array
        module procedure adios2_get_real_6d
        module procedure adios2_get_dp_6d
        module procedure adios2_get_complex_6d
        module procedure adios2_get_complex_dp_6d
        module procedure adios2_get_integer1_6d
        module procedure adios2_get_integer2_6d
        module procedure adios2_get_integer4_6d
        module procedure adios2_get_integer8_6d

        ! Single Value
        module procedure adios2_get_by_name_string
        module procedure adios2_get_by_name_real
        module procedure adios2_get_by_name_dp
        module procedure adios2_get_by_name_complex
        module procedure adios2_get_by_name_complex_dp
        module procedure adios2_get_by_name_integer1
        module procedure adios2_get_by_name_integer2
        module procedure adios2_get_by_name_integer4
        module procedure adios2_get_by_name_integer8

        ! 1D Array
        module procedure adios2_get_by_name_real_1d
        module procedure adios2_get_by_name_dp_1d
        module procedure adios2_get_by_name_complex_1d
        module procedure adios2_get_by_name_complex_dp_1d
        module procedure adios2_get_by_name_integer1_1d
        module procedure adios2_get_by_name_integer2_1d
        module procedure adios2_get_by_name_integer4_1d
        module procedure adios2_get_by_name_integer8_1d

        ! 2D Array
        module procedure adios2_get_by_name_real_2d
        module procedure adios2_get_by_name_dp_2d
        module procedure adios2_get_by_name_complex_2d
        module procedure adios2_get_by_name_complex_dp_2d
        module procedure adios2_get_by_name_integer1_2d
        module procedure adios2_get_by_name_integer2_2d
        module procedure adios2_get_by_name_integer4_2d
        module procedure adios2_get_by_name_integer8_2d

        ! 3D Array
        module procedure adios2_get_by_name_real_3d
        module procedure adios2_get_by_name_dp_3d
        module procedure adios2_get_by_name_complex_3d
        module procedure adios2_get_by_name_complex_dp_3d
        module procedure adios2_get_by_name_integer1_3d
        module procedure adios2_get_by_name_integer2_3d
        module procedure adios2_get_by_name_integer4_3d
        module procedure adios2_get_by_name_integer8_3d

        ! 4D Array
        module procedure adios2_get_by_name_real_4d
        module procedure adios2_get_by_name_dp_4d
        module procedure adios2_get_by_name_complex_4d
        module procedure adios2_get_by_name_complex_dp_4d
        module procedure adios2_get_by_name_integer1_4d
        module procedure adios2_get_by_name_integer2_4d
        module procedure adios2_get_by_name_integer4_4d
        module procedure adios2_get_by_name_integer8_4d

        ! 5D Array
        module procedure adios2_get_by_name_real_5d
        module procedure adios2_get_by_name_dp_5d
        module procedure adios2_get_by_name_complex_5d
        module procedure adios2_get_by_name_complex_dp_5d
        module procedure adios2_get_by_name_integer1_5d
        module procedure adios2_get_by_name_integer2_5d
        module procedure adios2_get_by_name_integer4_5d
        module procedure adios2_get_by_name_integer8_5d

        ! 6D Array
        module procedure adios2_get_by_name_real_6d
        module procedure adios2_get_by_name_dp_6d
        module procedure adios2_get_by_name_complex_6d
        module procedure adios2_get_by_name_complex_dp_6d
        module procedure adios2_get_by_name_integer1_6d
        module procedure adios2_get_by_name_integer2_6d
        module procedure adios2_get_by_name_integer4_6d
        module procedure adios2_get_by_name_integer8_6d

        ! Single Value
        module procedure adios2_get_deferred_string
        module procedure adios2_get_deferred_real
        module procedure adios2_get_deferred_dp
        module procedure adios2_get_deferred_complex
        module procedure adios2_get_deferred_complex_dp
        module procedure adios2_get_deferred_integer1
        module procedure adios2_get_deferred_integer2
        module procedure adios2_get_deferred_integer4
        module procedure adios2_get_deferred_integer8

        ! 1D Array
        module procedure adios2_get_deferred_real_1d
        module procedure adios2_get_deferred_dp_1d
        module procedure adios2_get_deferred_complex_1d
        module procedure adios2_get_deferred_complex_dp_1d
        module procedure adios2_get_deferred_integer1_1d
        module procedure adios2_get_deferred_integer2_1d
        module procedure adios2_get_deferred_integer4_1d
        module procedure adios2_get_deferred_integer8_1d

        ! 2D Array
        module procedure adios2_get_deferred_real_2d
        module procedure adios2_get_deferred_dp_2d
        module procedure adios2_get_deferred_complex_2d
        module procedure adios2_get_deferred_complex_dp_2d
        module procedure adios2_get_deferred_integer1_2d
        module procedure adios2_get_deferred_integer2_2d
        module procedure adios2_get_deferred_integer4_2d
        module procedure adios2_get_deferred_integer8_2d

        ! 3D Array
        module procedure adios2_get_deferred_real_3d
        module procedure adios2_get_deferred_dp_3d
        module procedure adios2_get_deferred_complex_3d
        module procedure adios2_get_deferred_complex_dp_3d
        module procedure adios2_get_deferred_integer1_3d
        module procedure adios2_get_deferred_integer2_3d
        module procedure adios2_get_deferred_integer4_3d
        module procedure adios2_get_deferred_integer8_3d

        ! 4D Array
        module procedure adios2_get_deferred_real_4d
        module procedure adios2_get_deferred_dp_4d
        module procedure adios2_get_deferred_complex_4d
        module procedure adios2_get_deferred_complex_dp_4d
        module procedure adios2_get_deferred_integer1_4d
        module procedure adios2_get_deferred_integer2_4d
        module procedure adios2_get_deferred_integer4_4d
        module procedure adios2_get_deferred_integer8_4d

        ! 5D Array
        module procedure adios2_get_deferred_real_5d
        module procedure adios2_get_deferred_dp_5d
        module procedure adios2_get_deferred_complex_5d
        module procedure adios2_get_deferred_complex_dp_5d
        module procedure adios2_get_deferred_integer1_5d
        module procedure adios2_get_deferred_integer2_5d
        module procedure adios2_get_deferred_integer4_5d
        module procedure adios2_get_deferred_integer8_5d

        ! 6D Array
        module procedure adios2_get_deferred_real_6d
        module procedure adios2_get_deferred_dp_6d
        module procedure adios2_get_deferred_complex_6d
        module procedure adios2_get_deferred_complex_dp_6d
        module procedure adios2_get_deferred_integer1_6d
        module procedure adios2_get_deferred_integer2_6d
        module procedure adios2_get_deferred_integer4_6d
        module procedure adios2_get_deferred_integer8_6d

        ! Deferred Signature
        ! Single Value
        module procedure adios2_get_deferred_by_name_string
        module procedure adios2_get_deferred_by_name_real
        module procedure adios2_get_deferred_by_name_dp
        module procedure adios2_get_deferred_by_name_complex
        module procedure adios2_get_deferred_by_name_complex_dp
        module procedure adios2_get_deferred_by_name_integer1
        module procedure adios2_get_deferred_by_name_integer2
        module procedure adios2_get_deferred_by_name_integer4
        module procedure adios2_get_deferred_by_name_integer8

        ! 1D Array
        module procedure adios2_get_deferred_by_name_real_1d
        module procedure adios2_get_deferred_by_name_dp_1d
        module procedure adios2_get_deferred_by_name_complex_1d
        module procedure adios2_get_deferred_by_name_complex_dp_1d
        module procedure adios2_get_deferred_by_name_integer1_1d
        module procedure adios2_get_deferred_by_name_integer2_1d
        module procedure adios2_get_deferred_by_name_integer4_1d
        module procedure adios2_get_deferred_by_name_integer8_1d

        ! 2D Array
        module procedure adios2_get_deferred_by_name_real_2d
        module procedure adios2_get_deferred_by_name_dp_2d
        module procedure adios2_get_deferred_by_name_complex_2d
        module procedure adios2_get_deferred_by_name_complex_dp_2d
        module procedure adios2_get_deferred_by_name_integer1_2d
        module procedure adios2_get_deferred_by_name_integer2_2d
        module procedure adios2_get_deferred_by_name_integer4_2d
        module procedure adios2_get_deferred_by_name_integer8_2d

        ! 3D Array
        module procedure adios2_get_deferred_by_name_real_3d
        module procedure adios2_get_deferred_by_name_dp_3d
        module procedure adios2_get_deferred_by_name_complex_3d
        module procedure adios2_get_deferred_by_name_complex_dp_3d
        module procedure adios2_get_deferred_by_name_integer1_3d
        module procedure adios2_get_deferred_by_name_integer2_3d
        module procedure adios2_get_deferred_by_name_integer4_3d
        module procedure adios2_get_deferred_by_name_integer8_3d

        ! 4D Array
        module procedure adios2_get_deferred_by_name_real_4d
        module procedure adios2_get_deferred_by_name_dp_4d
        module procedure adios2_get_deferred_by_name_complex_4d
        module procedure adios2_get_deferred_by_name_complex_dp_4d
        module procedure adios2_get_deferred_by_name_integer1_4d
        module procedure adios2_get_deferred_by_name_integer2_4d
        module procedure adios2_get_deferred_by_name_integer4_4d
        module procedure adios2_get_deferred_by_name_integer8_4d

        ! 5D Array
        module procedure adios2_get_deferred_by_name_real_5d
        module procedure adios2_get_deferred_by_name_dp_5d
        module procedure adios2_get_deferred_by_name_complex_5d
        module procedure adios2_get_deferred_by_name_complex_dp_5d
        module procedure adios2_get_deferred_by_name_integer1_5d
        module procedure adios2_get_deferred_by_name_integer2_5d
        module procedure adios2_get_deferred_by_name_integer4_5d
        module procedure adios2_get_deferred_by_name_integer8_5d

        ! 6D Array
        module procedure adios2_get_deferred_by_name_real_6d
        module procedure adios2_get_deferred_by_name_dp_6d
        module procedure adios2_get_deferred_by_name_complex_6d
        module procedure adios2_get_deferred_by_name_complex_dp_6d
        module procedure adios2_get_deferred_by_name_integer1_6d
        module procedure adios2_get_deferred_by_name_integer2_6d
        module procedure adios2_get_deferred_by_name_integer4_6d
        module procedure adios2_get_deferred_by_name_integer8_6d


    end interface
    external adios2_get_f2c
    external adios2_get_by_name_f2c

contains

    include 'contains/adios2_engine_get.f90'
    include 'contains/adios2_engine_get_by_name.f90'
    include 'contains/adios2_engine_get_deferred.f90'
    include 'contains/adios2_engine_get_deferred_by_name.f90'

end module
