program TestBPWriteAttributes
    use small_test_data
    use mpi
    use adios2
    implicit none

    type(adios2_adios) :: adios
    type(adios2_io) :: ioWrite, ioRead
    type(adios2_engine) :: bpWriter, bpReader
    type(adios2_attribute), dimension(14) :: attributes, attributes_in

    integer :: ierr, i
    character(len=:), allocatable :: attrName
    character(len=23):: iString_value
    integer(kind=1):: i8_value
    integer(kind=2):: i16_value
    integer(kind=4):: i32_value
    integer(kind=8):: i64_value
    real :: r32_value
    real(kind=8):: r64_value

    character(len=15), dimension(3):: iString_array
    integer(kind=1), dimension(3):: i8_array
    integer(kind=2), dimension(3):: i16_array
    integer(kind=4), dimension(3):: i32_array
    integer(kind=8), dimension(3):: i64_array
    real, dimension(3) :: r32_array
    real(kind=8), dimension(3):: r64_array

    type(adios2_namestruct) :: namestruct
    character(len=4096), dimension(:), allocatable :: attrnamelist
    

    ! Launch MPI
    INTEGER provided

    ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)

    ! Create adios handler passing the communicator and error flag
    call adios2_init(adios, MPI_COMM_WORLD, ierr)

    !!!!!!!!!!!!!!!!!!!!!!!!! WRITER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Declare an IO process configuration inside adios
    call adios2_declare_io(ioWrite, adios, "ioWrite", ierr)

    do i=1,14
       if( attributes(i)%valid .eqv. .true. ) then
          write(*,*) 'Invalid attribute default'
          stop 1
       end if
    end do

    ! single value
    call adios2_define_attribute(attributes(1), ioWrite, 'att_String', &
                                 'ADIOS2 String attribute', ierr)

    call adios2_define_attribute(attributes(2), ioWrite, 'att_i8', &
                                 data_I8(1), ierr)

    call adios2_define_attribute(attributes(3), ioWrite, 'att_i16', &
                                 data_I16(1), ierr)

    call adios2_define_attribute(attributes(4), ioWrite, 'att_i32', &
                                 data_I32(1), ierr)

    call adios2_define_attribute(attributes(5), ioWrite, 'att_i64', &
                                 data_I64(1), ierr)

    call adios2_define_attribute(attributes(6), ioWrite, 'att_r32', &
                                 data_R32(1), ierr)

    call adios2_define_attribute(attributes(7), ioWrite, 'att_r64', &
                                 data_R64(1), ierr)

    ! arrays
    call adios2_define_attribute(attributes(8), ioWrite, 'att_Strings_array', &
                                 data_Strings, 3, ierr)

    call adios2_define_attribute(attributes(9), ioWrite, 'att_i8_array', &
                                 data_I8, 3, ierr)

    call adios2_define_attribute(attributes(10), ioWrite, 'att_i16_array', &
                                 data_I16, 3, ierr)

    call adios2_define_attribute(attributes(11), ioWrite, 'att_i32_array', &
                                 data_I32, 3, ierr)

    call adios2_define_attribute(attributes(12), ioWrite, 'att_i64_array', &
                                 data_I64, 3, ierr)

    call adios2_define_attribute(attributes(13), ioWrite, 'att_r32_array', &
                                 data_R32, 3, ierr)

    call adios2_define_attribute(attributes(14), ioWrite, 'att_r64_array', &
                                 data_R64, 3, ierr)

    do i=1,14
       if( attributes(i)%valid .eqv. .false. ) then
          write(*,*) 'Invalid adios2_define_attribute'
          stop 1
       end if
    end do
     
    ! Testing adios2_attribute_name for just one case
    call adios2_attribute_name(attrName, attributes(1), ierr)
    if (attrName /= 'att_String') then
       write(*,*) 'Invalid adios2_attribute_name'
       stop 1
    end if
    deallocate(attrName)

    call adios2_open(bpWriter, ioWrite, "fattr_types.bp", adios2_mode_write, &
                     ierr)

    call adios2_close(bpWriter, ierr)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    !!!!!!!!!!!!!!!!!!!!!!!!! READER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call adios2_declare_io(ioRead, adios, "ioRead", ierr)

    call adios2_open(bpReader, ioRead, 'fattr_types.bp', adios2_mode_readRandomAccess, ierr)


   ! Test getting list of attribute names
    call adios2_available_attributes(ioRead, namestruct, ierr)
    if (ierr /= 0) then
       write(*,*) 'adios2_available_attributes returned with error'
       stop 1
    end if
    if (.not.namestruct%valid) then
       write(*,*) 'adios2_available_attributes returned invalid struct'
       stop 1
    end if
    write(*,*) 'Number of attributes = ', namestruct%count
    write(*,*) 'Max name length = ', namestruct%max_name_len
    if (namestruct%count /= 14) then
       write(*,*) 'adios2_available_attributes returned not the expected 14'
       stop 1
    end if

    allocate(attrnamelist(namestruct%count))

    call adios2_retrieve_names(namestruct, attrnamelist, ierr)
    if (ierr /= 0) then
       write(*,*) 'adios2_retrieve_names returned with error'
       stop 1
    end if
    do i=1,namestruct%count
         write(*,'("Attr[",i2,"] = ",a20)') i, attrnamelist(i)
    end do
    deallocate(attrnamelist)

    if (namestruct%f2c /= 0_8) then
       write(*,*) 'namestruct f2c pointer is not null after adios2_retrieve_names()'
       stop 1
    end if
    if (namestruct%valid) then
       write(*,*) 'namestruct is not invalidated after adios2_retrieve_names()'
       stop 1
    end if


    call adios2_inquire_attribute(attributes_in(1), ioRead, 'att_String', ierr)
    call adios2_inquire_attribute(attributes_in(2), ioRead, 'att_i8', ierr)
    call adios2_inquire_attribute(attributes_in(3), ioRead, 'att_i16', ierr)
    call adios2_inquire_attribute(attributes_in(4), ioRead, 'att_i32', ierr)
    call adios2_inquire_attribute(attributes_in(5), ioRead, 'att_i64', ierr)
    call adios2_inquire_attribute(attributes_in(6), ioRead, 'att_r32', ierr)
    call adios2_inquire_attribute(attributes_in(7), ioRead, 'att_r64', ierr)

    if(attributes_in(1)%valid .eqv. .false.) then
       write(*,*) 'attribute iString not found'
       stop 1
    end if
    if(attributes_in(1)%type /= adios2_type_string) then
       write(*,*) 'attribute iString wrong type'
       stop 1
    end if
    if(attributes_in(1)%length /= 1) then
       write(*,*) 'attribute iString length is not 1'
       stop 1
    end if
    if(attributes_in(1)%is_value .eqv. .false.) then
       write(*,*) 'attribute iString must be value'
       stop 1
    end if
    call adios2_attribute_data( iString_value, attributes_in(1), ierr)
    if( iString_value /=  'ADIOS2 String attribute' ) then
       write(*,*) 'attribute iString data error'
       stop 1
    end if

    if(attributes_in(2)%valid .eqv. .false.) then
       write(*,*) 'attribute i8 not found'
       stop 1
    end if
    if(attributes_in(2)%type /= adios2_type_integer1) then
       write(*,*) 'attribute i8 wrong type'
       stop 1
    end if
    if(attributes_in(2)%length /= 1) then
       write(*,*) 'attribute i8 length is not 1'
       stop 1
    end if
    if(attributes_in(2)%is_value .eqv. .false.) then
       write(*,*) 'attribute i8 must be value'
       stop 1
    end if
    call adios2_attribute_data( i8_value, attributes_in(2), ierr)
    if( i8_value /=  data_I8(1) ) then
       write(*,*) 'attribute i8 data error'
       stop 1
    end if

    if(attributes_in(3)%valid .eqv. .false.) then
       write(*,*) 'attribute i16 not found'
       stop 1
    end if
    if(attributes_in(3)%type /= adios2_type_integer2) then
       write(*,*) 'attribute i16 wrong type'
       stop 1
    end if
    if(attributes_in(3)%length /= 1) then
       write(*,*) 'attribute i16 length is not 1'
       stop 1
    end if
    if(attributes_in(3)%is_value .eqv. .false.) then
       write(*,*) 'attribute i16 must be value'
       stop 1
    end if
    call adios2_attribute_data( i16_value, attributes_in(3), ierr)
    if( i16_value /=  data_I16(1) ) then
       write(*,*) 'attribute i16 data error'
       stop 1
    end if

    if(attributes_in(4)%valid .eqv. .false.) then
       write(*,*) 'attribute i32 not found'
       stop 1
    end if
    if(attributes_in(4)%type /= adios2_type_integer4) then
       write(*,*) 'attribute i32 wrong type'
       stop 1
    end if
    if(attributes_in(4)%length /= 1) then
       write(*,*) 'attribute i32 length is not 1'
       stop 1
    end if
    if(attributes_in(4)%is_value .eqv. .false.) then
       write(*,*) 'attribute i32 must be value'
       stop 1
    end if
    call adios2_attribute_data( i32_value, attributes_in(4), ierr)
    if( i32_value /=  data_I32(1) ) then
       write(*,*) 'attribute i32 data error'
       stop 1
    end if

    if(attributes_in(5)%valid .eqv. .false.) then
       write(*,*) 'attribute i64 not found'
       stop 1
    end if
    if(attributes_in(5)%type /= adios2_type_integer8) then
       write(*,*) 'attribute i64 wrong type'
       stop 1
    end if
    if(attributes_in(5)%length /= 1) then
       write(*,*) 'attribute i64 length is not 1'
       stop 1
    end if
    if(attributes_in(5)%is_value .eqv. .false.) then
       write(*,*) 'attribute i64 must be value'
       stop 1
    end if
    call adios2_attribute_data( i64_value, attributes_in(5), ierr)
    if( i64_value /=  data_I64(1) ) then
       write(*,*) 'attribute i64 data error'
       stop 1
    end if

    if(attributes_in(6)%valid .eqv. .false.) then
       write(*,*) 'attribute r32 not found'
       stop 1
    end if
    if(attributes_in(6)%type /= adios2_type_real) then
       write(*,*) 'attribute r32 wrong type'
       stop 1
    end if
    if(attributes_in(6)%length /= 1) then
       write(*,*) 'attribute r32 length is not 1'
       stop 1
    end if
    if(attributes_in(6)%is_value .eqv. .false.) then
       write(*,*) 'attribute r32 must be value'
       stop 1
    end if
    call adios2_attribute_data( r32_value, attributes_in(6), ierr)
    if( r32_value /=  data_R32(1) ) then
       write(*,*) 'attribute r32 data error'
       stop 1
    end if

    if(attributes_in(7)%valid .eqv. .false.) then
       write(*,*) 'attribute r64 not found'
       stop 1
    end if
    if(attributes_in(7)%type /= adios2_type_dp) then
       write(*,*) 'attribute r64 wrong type'
       stop 1
    end if
    if(attributes_in(7)%length /= 1) then
       write(*,*) 'attribute r64 length is not 1'
       stop 1
    end if
    if(attributes_in(7)%is_value .eqv. .false.) then
       write(*,*) 'attribute r64 must be value'
       stop 1
    end if
    call adios2_attribute_data( r64_value, attributes_in(7), ierr)
    if( r64_value /=  data_R64(1) ) then
       write(*,*) 'attribute r64 data error'
       stop 1
    end if

    ! Array
    call adios2_inquire_attribute(attributes_in(8), ioRead, 'att_Strings_array', ierr)
    call adios2_inquire_attribute(attributes_in(9), ioRead, 'att_i8_array', ierr)
    call adios2_inquire_attribute(attributes_in(10), ioRead, 'att_i16_array', ierr)
    call adios2_inquire_attribute(attributes_in(11), ioRead, 'att_i32_array', ierr)
    call adios2_inquire_attribute(attributes_in(12), ioRead, 'att_i64_array', ierr)
    call adios2_inquire_attribute(attributes_in(13), ioRead, 'att_r32_array', ierr)
    call adios2_inquire_attribute(attributes_in(14), ioRead, 'att_r64_array', ierr)

    if(attributes_in(8)%valid .eqv. .false.) then
       write(*,*) 'attribute string array not found'
       stop 1
    end if
    if(attributes_in(8)%type /= adios2_type_string) then
       write(*,*) 'attribute string array wrong type'
       stop 1
    end if
    if(attributes_in(8)%length /= 3) then
       write(*,*) 'attribute string array length is not 3'
       stop 1
    end if
    if(attributes_in(8)%is_value .eqv. .true.) then
       write(*,*) 'attribute string array must be array'
       stop 1
    end if
    call adios2_attribute_data( iString_array, attributes_in(8), ierr)
    do i=1,3
       if( iString_array(i) /= data_Strings(i) ) then
          write(*,*) 'attribute string array data error'
          stop 1
       end if
    end do

    if(attributes_in(9)%valid .eqv. .false.) then
       write(*,*) 'attribute i8 array not found'
       stop 1
    end if
    if(attributes_in(9)%type /= adios2_type_integer1) then
       write(*,*) 'attribute i8 array wrong type'
       stop 1
    end if
    if(attributes_in(9)%length /= 3) then
       write(*,*) 'attribute i8 array length is not 3'
       stop 1
    end if
    if(attributes_in(9)%is_value .eqv. .true.) then
       write(*,*) 'attribute i8 array must be array'
       stop 1
    end if
    call adios2_attribute_data( i8_array, attributes_in(9), ierr)
    do i=1,3
       if( i8_array(i) /=  data_I8(i) ) then
          write(*,*) 'attribute i8 array data error'
          stop 1
       end if
    end do

    if(attributes_in(10)%valid .eqv. .false.) then
       write(*,*) 'attribute i16 array not found'
       stop 1
    end if
    if(attributes_in(10)%type /= adios2_type_integer2) then
       write(*,*) 'attribute i16 array wrong type'
       stop 1
    end if
    if(attributes_in(10)%length /= 3) then
       write(*,*) 'attribute i16 array length is not 3'
       stop 1
    end if
    if(attributes_in(10)%is_value .eqv. .true.) then
       write(*,*) 'attribute i16 array must be array'
       stop 1
    end if
    call adios2_attribute_data( i16_array, attributes_in(10), ierr)
    do i=1,3
       if( i16_array(i) /=  data_I16(i) ) then
          write(*,*) 'attribute i16 array data error'
          stop 1
       end if
    end do

    if(attributes_in(11)%valid .eqv. .false.) then
       write(*,*) 'attribute i32 array not found'
       stop 1
    end if
    if(attributes_in(11)%type /= adios2_type_integer4) then
       write(*,*) 'attribute i32 array wrong type'
       stop 1
    end if
    if(attributes_in(11)%length /= 3) then
       write(*,*) 'attribute i32 array length is not 3'
       stop 1
    end if
    if(attributes_in(11)%is_value .eqv. .true.) then
       write(*,*) 'attribute i32 array must be array'
       stop 1
    end if
    call adios2_attribute_data( i32_array, attributes_in(11), ierr)
    do i=1,3
       if( i32_array(i) /=  data_I32(i) ) then
          write(*,*) 'attribute i32 array data error'
          stop 1
       end if
    end do

    if(attributes_in(12)%valid .eqv. .false.) then
       write(*,*) 'attribute i64 array not found'
       stop 1
    end if
    if(attributes_in(12)%type /= adios2_type_integer8) then
       write(*,*) 'attribute i64 array wrong type'
       stop 1
    end if
    if(attributes_in(12)%length /= 3) then
       write(*,*) 'attribute i64 array length is not 3'
       stop 1
    end if
    if(attributes_in(12)%is_value .eqv. .true.) then
       write(*,*) 'attribute i64 array must be array'
       stop 1
    end if
    call adios2_attribute_data( i64_array, attributes_in(12), ierr)
    do i=1,3
       if( i64_array(i) /=  data_I64(i) ) then
          write(*,*) 'attribute i64 array data error'
          stop 1
       end if
    end do

    if(attributes_in(13)%valid .eqv. .false.) then
       write(*,*) 'attribute r32 array not found'
       stop 1
    end if
    if(attributes_in(13)%type /= adios2_type_real) then
       write(*,*) 'attribute r32 array wrong type'
       stop 1
    end if
    if(attributes_in(13)%length /= 3) then
       write(*,*) 'attribute r32 array length is not 3'
       stop 1
    end if
    if(attributes_in(13)%is_value .eqv. .true.) then
       write(*,*) 'attribute r32 array must be array'
       stop 1
    end if
    call adios2_attribute_data( r32_array, attributes_in(13), ierr)
    do i=1,3
       if( r32_array(i) /=  data_R32(i) ) then
          write(*,*) 'attribute r32 array data error'
          stop 1
       end if
    end do

    if(attributes_in(14)%valid .eqv. .false.) then
       write(*,*) 'attribute r64 array not found'
       stop 1
    end if
    if(attributes_in(14)%type /= adios2_type_dp) then
       write(*,*) 'attribute r64 array wrong type'
       stop 1
    end if
    if(attributes_in(14)%length /= 3) then
       write(*,*) 'attribute r64 array length is not 3'
       stop 1
    end if
    if(attributes_in(14)%is_value .eqv. .true.) then
       write(*,*) 'attribute r64 array must be array'
       stop 1
    end if
    call adios2_attribute_data( r64_array, attributes_in(14), ierr)
    do i=1,3
       if( r64_array(i) /=  data_R64(i) ) then
          write(*,*) 'attribute r64 array data error'
          stop 1
       end if
    end do

    call adios2_close(bpReader, ierr)



    call adios2_finalize(adios, ierr)

    call MPI_Finalize(ierr)

end program TestBPWriteAttributes
