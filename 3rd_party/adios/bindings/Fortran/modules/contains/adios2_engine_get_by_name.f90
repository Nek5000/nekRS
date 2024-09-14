!
! Distributed under the OSI-approved Apache License, Version 2.0.  See
!  accompanying file Copyright.txt for details.
!
!  adios2_engine_get_by_name.f90 : implementation of adios2_get_by_name
!  subroutines
!
!   Created on: Feb 21, 2018
!       Author: William F Godoy godoywf@ornl.gov
!

! Single values
subroutine adios2_get_by_name_string(engine, name, data, launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    character*(*), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_real(engine, name, data, launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real,  intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_dp(engine, name, data, &
                                    launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real(kind=8),  intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex(engine, name, data, &
                                         launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex,  intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_dp(engine, name, data, &
                                            launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex(kind=8),  intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer1(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=1),  intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer2(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=2),  intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer4(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=4),  intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer8(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=8),  intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

! 1D arrays
subroutine adios2_get_by_name_real_1d(engine, name, data, &
                                      launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real, dimension(:), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_dp_1d(engine, name, data, &
                                    launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real(kind=8), dimension(:), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_1d(engine, name, data, &
                                         launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex, dimension(:), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_dp_1d(engine, name, data, &
                                            launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex(kind=8), dimension(:), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer1_1d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=1), dimension(:), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer2_1d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=2), dimension(:), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer4_1d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=4), dimension(:), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer8_1d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=8), dimension(:), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

! 2D arrays
subroutine adios2_get_by_name_real_2d(engine, name, data, &
                                      launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real, dimension(:, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_dp_2d(engine, name, data, &
                                    launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real(kind=8), dimension(:, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_2d(engine, name, data, &
                                         launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex, dimension(:, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_dp_2d(engine, name, data, &
                                            launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex(kind=8), dimension(:, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer1_2d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=1), dimension(:, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer2_2d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=2), dimension(:, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer4_2d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=4), dimension(:, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer8_2d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=8), dimension(:, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

! 3D arrays
subroutine adios2_get_by_name_real_3d(engine, name, data, &
                                      launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real, dimension(:, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_dp_3d(engine, name, data, &
                                    launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real(kind=8), dimension(:, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_3d(engine, name, data, &
                                         launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex, dimension(:, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_dp_3d(engine, name, data, &
                                            launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex(kind=8), dimension(:, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer1_3d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=1), dimension(:, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer2_3d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=2), dimension(:, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer4_3d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=4), dimension(:, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer8_3d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=8), dimension(:, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

! 4D arrays
subroutine adios2_get_by_name_real_4d(engine, name, data, &
                                      launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real, dimension(:, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_dp_4d(engine, name, data, &
                                    launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real(kind=8), dimension(:, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_4d(engine, name, data, &
                                         launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex, dimension(:, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_dp_4d(engine, name, data, &
                                            launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex(kind=8), dimension(:, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer1_4d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=1), dimension(:, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer2_4d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=2), dimension(:, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer4_4d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=4), dimension(:, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer8_4d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=8), dimension(:, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

! 5D arrays
subroutine adios2_get_by_name_real_5d(engine, name, data, &
                                      launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real, dimension(:, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_dp_5d(engine, name, data, &
                                    launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real(kind=8), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_5d(engine, name, data, &
                                         launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex, dimension(:, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_dp_5d(engine, name, data, &
                                            launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex(kind=8), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer1_5d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=1), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer2_5d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=2), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer4_5d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=4), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer8_5d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=8), dimension(:, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

! 6D arrays
subroutine adios2_get_by_name_real_6d(engine, name, data, &
                                      launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real, dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_dp_6d(engine, name, data, &
                                    launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    real(kind=8), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_6d(engine, name, data, &
                                         launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex, dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_complex_dp_6d(engine, name, data, &
                                            launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    complex(kind=8), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer1_6d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=1), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer2_6d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=2), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer4_6d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=4), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

subroutine adios2_get_by_name_integer8_6d(engine, name, data, &
                                          launch, ierr)
    type(adios2_engine), intent(in):: engine
    character*(*), intent(in) :: name
    integer(kind=8), dimension(:, :, :, :, :, :), intent(out):: data
    integer, intent(in):: launch
    integer, intent(out):: ierr

    if(trim(engine%type) == "NULL") return
    call adios2_get_by_name_f2c(engine%f2c, &
                                TRIM(ADJUSTL(name))//char(0), &
                                data, launch, ierr)

end subroutine

