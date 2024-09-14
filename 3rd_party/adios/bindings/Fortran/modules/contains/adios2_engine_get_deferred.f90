!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_engine_get_deferred.f90 : implementation of adios2_get_deferred
!  subroutines
!
!   Created on: Feb 21, 2018
!       Author: William F Godoy godoywf@ornl.gov
!

! Single data
subroutine adios2_get_deferred_string(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    character*(*), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_string, &
                                    'get string', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, adios2_mode_sync, &
                            ierr)
    end if

end subroutine

subroutine adios2_get_deferred_real(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real, intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_real, &
                                    'get real', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_dp(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real(kind=8), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_dp, &
                                    'get dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex, intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex, &
                                    'get complex', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_dp(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex(kind=8), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex_dp, &
                                    'get complex_dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer1(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=1), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer1, &
                                    'get integer1', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer2(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=2), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer2, &
                                    'get integer2', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer4(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=4), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer4, &
                                    'get integer4', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer8(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=8), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer8, &
                                    'get integer8', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

! 1D arrays
subroutine adios2_get_deferred_real_1d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real, dimension(:), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_real, &
                                    'get real', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_dp_1d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real(kind=8), dimension(:), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_dp, &
                                    'get dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_1d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex, dimension(:), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex, &
                                    'get complex', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_dp_1d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex(kind=8), dimension(:), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex_dp, &
                                    'get complex_dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer1_1d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=1), dimension(:), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer1, &
                                    'get integer1', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer2_1d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=2), dimension(:), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer2, &
                                    'get integer2', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer4_1d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=4), dimension(:), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer4, &
                                    'get integer4', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer8_1d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=8), dimension(:), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer8, &
                                    'get integer8', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

! 2D arrays
subroutine adios2_get_deferred_real_2d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real, dimension(:, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_real, &
                                    'get real', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_dp_2d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real(kind=8), dimension(:, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_dp, &
                                    'get dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_2d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex, dimension(:, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex, &
                                    'get complex', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_dp_2d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex(kind=8), dimension(:, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex_dp, &
                                    'get complex_dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer1_2d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=1), dimension(:, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer1, &
                                    'get integer1', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer2_2d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=2), dimension(:, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer2, &
                                    'get integer2', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer4_2d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=4), dimension(:, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer4, &
                                    'get integer4', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer8_2d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=8), dimension(:, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer8, &
                                    'get integer8', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

! 3D arrays
subroutine adios2_get_deferred_real_3d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real, dimension(:, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_real, &
                                    'get real', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_dp_3d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real(kind=8), dimension(:, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_dp, &
                                    'get dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_3d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex, dimension(:, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex, &
                                    'get complex', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_dp_3d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex(kind=8), dimension(:, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex_dp, &
                                    'get complex_dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer1_3d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=1), dimension(:, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer1, &
                                    'get integer1', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer2_3d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=2), dimension(:, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer2, &
                                    'get integer2', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer4_3d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=4), dimension(:, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer4, &
                                    'get integer4', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer8_3d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=8), dimension(:, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer8, &
                                    'get integer8', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

! 4D arrays
subroutine adios2_get_deferred_real_4d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real, dimension(:, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_real, &
                                    'get real', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_dp_4d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real(kind=8), dimension(:, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_dp, &
                                    'get dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_4d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex, dimension(:, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex, &
                                    'get complex', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_dp_4d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex(kind=8), dimension(:, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex_dp, &
                                    'get complex_dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer1_4d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=1), dimension(:, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer1, &
                                    'get integer1', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer2_4d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=2), dimension(:, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer2, &
                                    'get integer2', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer4_4d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=4), dimension(:, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer4, &
                                    'get integer4', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer8_4d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=8), dimension(:, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer8, &
                                    'get integer8', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

! 5D arrays
subroutine adios2_get_deferred_real_5d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real, dimension(:, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_real, &
                                    'get real', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_dp_5d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real(kind=8), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_dp, &
                                    'get dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_5d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex, dimension(:, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex, &
                                    'get complex', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_dp_5d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex(kind=8), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex_dp, &
                                    'get complex_dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer1_5d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=1), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer1, &
                                    'get integer1', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer2_5d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=2), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer2, &
                                    'get integer2', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer4_5d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=4), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer4, &
                                    'get integer4', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer8_5d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=8), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer8, &
                                    'get integer8', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

! 6D arrays
subroutine adios2_get_deferred_real_6d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real, dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_real, &
                                    'get real', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_dp_6d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    real(kind=8), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_dp, &
                                    'get dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_6d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex, dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex, &
                                    'get complex', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_complex_dp_6d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    complex(kind=8), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_complex_dp, &
                                    'get complex_dp', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer1_6d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=1), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer1, &
                                    'get integer1', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer2_6d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=2), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer2, &
                                    'get integer2', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer4_6d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=4), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer4, &
                                    'get integer4', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine

subroutine adios2_get_deferred_integer8_6d(engine, variable, data, ierr)
    type(adios2_engine), intent(in):: engine
    type(adios2_variable), intent(in):: variable
    integer(kind=8), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_variable_check_type(variable, adios2_type_integer8, &
                                    'get integer8', ierr)
    if (ierr == 0) then
        call adios2_get_f2c(engine%f2c, variable%f2c, data, &
                            adios2_mode_deferred, ierr)
    end if

end subroutine
