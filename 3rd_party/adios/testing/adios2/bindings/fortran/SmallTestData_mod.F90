! Distributed under the OSI-approved Apache License, Version 2.0.  See
! accompanying file Copyright.txt for details.
!
! SmallTestData_mod.f90 : small Fortran 90 arrays data for tests
!
!  Created on: Aug 9, 2017
!      Author: William F Godoy godoywf@ornl.gov
!

module small_test_data
    implicit none

    character(len=15), parameter, dimension(3) :: data_Strings = &
                                                  (/'Attribute oneXX', 'Attribute twoXX', 'Attribute three'/)

    integer(kind=1), parameter, dimension(10) :: data_I8 = &
                                                 INT((/0, 1, -2, 3, -4, 5, -6, 7, -8, 9/), 1)

    integer(kind=2), parameter, dimension(10) :: data_I16 = &
                                                 INT((/512, 513, -510, 515, -508, 517, -506, 519, -504, 521/), 2)

    integer(kind=4), parameter, dimension(10) :: data_I32 = &
                                                 (/131072, 131073, -131070, 131075, -131068, 131077, -131066, 131079, &
                                                   -131064, 131081/)

    integer(kind=8), parameter, dimension(10) :: data_I64 = &
                                                 (/8589934592_8, 8589934593_8, -8589934590_8, 8589934595_8, &
                                                   -8589934588_8, 8589934597_8, -8589934586_8, 8589934599_8, &
                                                   -8589934584_8, 8589934601_8/)

    real(kind=4), parameter, dimension(10) :: data_R32 = &
                                              (/0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1/)

    real(kind=8), parameter, dimension(10) :: data_R64 = &
                                              (/10.2, 11.2, 12.2, 13.2, 14.2, 15.2, 16.2, 17.2, 18.2, 19.2/)

end module
