program TestBPWriteReadHeatMapZfp2DRemove
  use mpi
  use adios2

  implicit none

  type(adios2_adios) :: adios
  type(adios2_io) :: ioPut, ioGet
  type(adios2_engine) :: bpWriter, bpReader
  type(adios2_variable), dimension(2) :: var_temperatures, var_temperaturesIn
  type(adios2_operator) :: zfp_operator
  integer:: operation_id
  integer(kind=4) :: step_status, i

  real(kind=4), dimension(:, :), allocatable :: temperatures_r4, &
                                                sel_temperatures_r4

  real(kind=8), dimension(:, :), allocatable :: temperatures_r8, &
                                                sel_temperatures_r8

  integer(kind=8), dimension(2) :: ishape, istart, icount
  integer(kind=8), dimension(2) :: sel_start, sel_count
  integer :: ierr, irank, isize
  integer :: in1, in2

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)

  in1 = 10
  in2 = 10

  icount = (/in1, in2/)
  istart = (/0, in2*irank/)
  ishape = (/in1, in2*isize/)

  allocate (temperatures_r4(in1, in2))
  allocate (temperatures_r8(in1, in2))

  temperatures_r4 = 1.0
  temperatures_r8 = 1.0_8

  ! Start adios2 Writer
  call adios2_init(adios, MPI_COMM_WORLD, ierr)
  
  call adios2_define_operator(zfp_operator, adios, 'CompressorZfp', 'zfp', ierr)
  call adios2_declare_io(ioPut, adios, 'HeatMapWrite', ierr)
  
  call adios2_define_variable(var_temperatures(1), ioPut, &
                              'temperatures_r4', adios2_type_real, &
                              2, ishape, istart, icount, &
                              adios2_constant_dims, ierr)

  call adios2_define_variable(var_temperatures(2), ioPut, &
                              'temperatures_r8', adios2_type_dp, &
                              2, ishape, istart, icount, &
                              adios2_constant_dims, ierr)

  call adios2_open(bpWriter, ioPut, 'HeatMapZfp2DRemove_f.bp', adios2_mode_write, &
                   ierr)

  do i=0,4

    call adios2_begin_step(bpWriter, adios2_step_mode_append, -1.0, &
                           step_status, ierr)
    
    if( mod(i,2) == 0 ) then
      call adios2_add_operation(operation_id, var_temperatures(1), &
                            zfp_operator, 'rate', '8', ierr)
      if( operation_id /= 0 ) then
         write(*,*) 'operation_id not added for real type'
         stop 1
      end if
      
      call adios2_add_operation(operation_id, var_temperatures(2), &
                            zfp_operator, 'rate', '8', ierr)
      if( operation_id /= 0 ) then
         write(*,*) 'operation_id not added for dp type'
         stop 1
      end if
                    
    else
      call adios2_remove_operations(var_temperatures(1), ierr)
      call adios2_remove_operations(var_temperatures(2), ierr)
    end if
    
    call adios2_put(bpWriter, var_temperatures(1), temperatures_r4, ierr)
    call adios2_put(bpWriter, var_temperatures(2), temperatures_r8, ierr)

    call adios2_end_step(bpWriter, ierr)
  end do

  call adios2_close(bpWriter, ierr)

  if (allocated(temperatures_r4)) deallocate (temperatures_r4)
  if (allocated(temperatures_r8)) deallocate (temperatures_r8)

  ! Start adios2 Reader in rank 0
  if (irank == 0) then

    call adios2_declare_io(ioGet, adios, 'HeatMapRead', ierr)

    call adios2_open(bpReader, ioGet, 'HeatMapZfp2DRemove_f.bp', &
                     adios2_mode_read, MPI_COMM_SELF, ierr)

    sel_start = (/0, 0/)
    sel_count = (/ishape(1), ishape(2)/)

    allocate (sel_temperatures_r4(sel_count(1), sel_count(2)))
    allocate (sel_temperatures_r8(sel_count(1), sel_count(2)))

    sel_temperatures_r4 = 0.0_4
    sel_temperatures_r8 = 0.0_8
    
    do i=0,4
      
      call adios2_begin_step(bpReader, adios2_step_mode_read, -1.0, &
                             step_status, ierr)  

      call adios2_inquire_variable(var_temperaturesIn(1), ioGet, &
                                 'temperatures_r4', ierr)
      call adios2_inquire_variable(var_temperaturesIn(2), ioGet, &
                                 'temperatures_r8', ierr)

      call adios2_set_selection(var_temperaturesIn(1), 2, sel_start, sel_count, &
                              ierr)
      call adios2_set_selection(var_temperaturesIn(2), 2, sel_start, sel_count, &
                              ierr)
  
      call adios2_get(bpReader, var_temperaturesIn(1), sel_temperatures_r4, ierr)
      call adios2_get(bpReader, var_temperaturesIn(2), sel_temperatures_r8, ierr)
      call adios2_end_step(bpReader, ierr)

      if (sum(sel_temperatures_r4) /= 100*isize) then
         write(*,*) 'Test failed real*4'
         stop 1
      end if
      if (sum(sel_temperatures_r8) /= 100*isize) then
         write(*,*) 'Test failed real*8'
         stop 1
      end if
    
    end do  

    call adios2_close(bpReader, ierr)

    if (allocated(sel_temperatures_r4)) deallocate (sel_temperatures_r4)
    if (allocated(sel_temperatures_r8)) deallocate (sel_temperatures_r8)

  end if

  call adios2_finalize(adios, ierr)
  call MPI_Finalize(ierr)

end program TestBPWriteReadHeatMapZfp2DRemove
