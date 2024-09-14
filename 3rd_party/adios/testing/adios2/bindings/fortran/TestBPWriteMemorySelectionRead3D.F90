program TestBPWriteMemorySelectionRead3D
  use mpi
  use adios2

  implicit none

  type(adios2_adios) :: adios
  type(adios2_io) :: ioPut, ioGet
  type(adios2_engine) :: bpWriter, bpReader
  type(adios2_variable), dimension(6) :: vars, vars_in

  integer(kind=1), dimension(:, :, :), allocatable :: data_i1, in_data_i1
  integer(kind=2), dimension(:, :, :), allocatable :: data_i2, in_data_i2
  integer(kind=4), dimension(:, :, :), allocatable :: data_i4, in_data_i4
  integer(kind=8), dimension(:, :, :), allocatable :: data_i8, in_data_i8
  real(kind=4), dimension(:, :, :), allocatable :: data_r4, in_data_r4
  real(kind=8), dimension(:, :, :), allocatable :: data_r8, in_data_r8

  integer(kind=1) :: min_i1, max_i1
  integer(kind=2) :: min_i2, max_i2
  integer(kind=4) :: min_i4, max_i4
  integer(kind=8) :: min_i8, max_i8
  real(kind=4)    :: min_r4, max_r4
  real(kind=8)    :: min_r8, max_r8

  integer(kind=8), dimension(3) :: ishape, istart, icount
  integer(kind=8), dimension(3) :: memory_start, memory_count
  integer :: ierr, irank, isize, nx, ny, nz, nsteps, s, step_status
  integer(kind=8) :: ghost_x, ghost_y, ghost_z, i, j, k, current_step
  integer(kind=8) :: expected_min, expected_max

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)

  nsteps = 3

  nx = 10 ! fastest
  ny = 20
  nz = 2

  ghost_x = 2
  ghost_y = 4
  ghost_z = 1

  allocate (data_i1(nx + 2* ghost_x, ny + 2* ghost_y, nz + 2*ghost_z))
  allocate (data_i2(nx + 2* ghost_x, ny + 2* ghost_y, nz + 2*ghost_z))
  allocate (data_i4(nx + 2* ghost_x, ny + 2* ghost_y, nz + 2*ghost_z))
  allocate (data_i8(nx + 2* ghost_x, ny + 2* ghost_y, nz + 2*ghost_z))
  allocate (data_r4(nx + 2* ghost_x, ny + 2* ghost_y, nz + 2*ghost_z))
  allocate (data_r8(nx + 2* ghost_x, ny + 2* ghost_y, nz + 2*ghost_z))

  data_i1 = -1
  data_i2 = -1
  data_i4 = -1
  data_i8 = -1_8
  data_r4 = -1.0
  data_r8 = -1.0_8

  call adios2_init(adios, MPI_COMM_WORLD, ierr)

  ! Writer
  call adios2_declare_io(ioPut, adios, 'MemSelWriter', ierr)

  ishape = (/nx, ny, isize*nz/)
  istart = (/0,   0, irank*nz/)
  icount = (/nx, ny,       nz/)

  memory_start = (/ ghost_x, ghost_y, ghost_z /)
  memory_count = (/ nx + 2* ghost_x, ny + 2* ghost_y, nz + 2* ghost_z /)

  call adios2_define_variable(vars(1), ioPut, 'var_i1', adios2_type_integer1, &
                              3, ishape, istart, icount, &
                              adios2_constant_dims, ierr)

  call adios2_define_variable(vars(2), ioPut, 'var_i2', adios2_type_integer2, &
                              3, ishape, istart, icount, &
                              adios2_constant_dims, ierr)

  call adios2_define_variable(vars(3), ioPut, 'var_i4', adios2_type_integer4, &
                              3, ishape, istart, icount, &
                              adios2_constant_dims, ierr)

  call adios2_define_variable(vars(4), ioPut, 'var_i8', adios2_type_integer8, &
                              3, ishape, istart, icount, &
                              adios2_constant_dims, ierr)

  call adios2_define_variable(vars(5), ioPut, 'var_r4', adios2_type_real, &
                              3, ishape, istart, icount, &
                              adios2_constant_dims, ierr)

  call adios2_define_variable(vars(6), ioPut, 'var_r8', adios2_type_dp, &
                              3, ishape, istart, icount, &
                              adios2_constant_dims, ierr)

  call adios2_set_memory_selection(vars(1), 3, memory_start, memory_count, ierr)
  call adios2_set_memory_selection(vars(2), 3, memory_start, memory_count, ierr)
  call adios2_set_memory_selection(vars(3), 3, memory_start, memory_count, ierr)
  call adios2_set_memory_selection(vars(4), 3, memory_start, memory_count, ierr)
  call adios2_set_memory_selection(vars(5), 3, memory_start, memory_count, ierr)
  call adios2_set_memory_selection(vars(6), 3, memory_start, memory_count, ierr)

  call adios2_open(bpWriter, ioPut, 'MemSel3D_f.bp', adios2_mode_write, ierr)

  do s=0,nsteps-1
    data_i1(ghost_x+1:nx+ghost_x+1,ghost_y+1:ny+ghost_y+1,ghost_z+1:nz+ghost_z+1) = INT(s, 1)
    data_i2(ghost_x+1:nx+ghost_x+1,ghost_y+1:ny+ghost_y+1,ghost_z+1:nz+ghost_z+1) = INT(s, 2)
    data_i4(ghost_x+1:nx+ghost_x+1,ghost_y+1:ny+ghost_y+1,ghost_z+1:nz+ghost_z+1) = INT(s, 4)
    data_i8(ghost_x+1:nx+ghost_x+1,ghost_y+1:ny+ghost_y+1,ghost_z+1:nz+ghost_z+1) = INT(s, 8)
    data_r4(ghost_x+1:nx+ghost_x+1,ghost_y+1:ny+ghost_y+1,ghost_z+1:nz+ghost_z+1) = INT(s, 4)
    data_r8(ghost_x+1:nx+ghost_x+1,ghost_y+1:ny+ghost_y+1,ghost_z+1:nz+ghost_z+1) = INT(s, 8)

    call adios2_begin_step(bpWriter, ierr)
    call adios2_put(bpWriter, vars(1), data_i1, ierr)
    call adios2_put(bpWriter, vars(2), data_i2, ierr)
    call adios2_put(bpWriter, vars(3), data_i4, ierr)
    call adios2_put(bpWriter, vars(4), data_i8, ierr)
    call adios2_put(bpWriter, vars(5), data_r4, ierr)
    call adios2_put(bpWriter, vars(6), data_r8, ierr)
    call adios2_end_step(bpWriter, ierr)

  end do

  call adios2_close(bpWriter, ierr)

  deallocate (data_i1)
  deallocate (data_i2)
  deallocate (data_i4)
  deallocate (data_i8)
  deallocate (data_r4)
  deallocate (data_r8)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! reader goes here
  call adios2_declare_io(ioGet, adios, 'MemSelReader', ierr)

  call adios2_open(bpReader, ioGet, 'MemSel3D_f.bp', adios2_mode_read, ierr)

  allocate (in_data_i1(nx, ny, isize*nz))
  allocate (in_data_i2(nx, ny, isize*nz))
  allocate (in_data_i4(nx, ny, isize*nz))
  allocate (in_data_i8(nx, ny, isize*nz))
  allocate (in_data_r4(nx, ny, isize*nz))
  allocate (in_data_r8(nx, ny, isize*nz))

  do
    call adios2_begin_step(bpReader, adios2_step_mode_read, -1., &
                           step_status, ierr)

    if(step_status == adios2_step_status_end_of_stream) exit

    call adios2_current_step(current_step, bpReader, ierr)

    expected_min = current_step
    expected_max = current_step

    call adios2_inquire_variable(vars_in(1), ioGet, 'var_i1', ierr)
    if( vars_in(1)%valid .eqv. .false. ) then
       write(*,*) 'i1 invalid error'
       stop 1
    end if
    if( vars_in(1)%name /= 'var_i1' ) then
       write(*,*) 'i1 name error'
       stop 1
    end if
    if( vars_in(1)%type /= adios2_type_integer1 ) then
       write(*,*) 'i1 type error'
       stop 1
    end if
    if( vars_in(1)%ndims /= 3 ) then
       write(*,*) 'i1 dims error'
       stop 1
    end if
    call adios2_variable_min(min_i1, vars_in(1), ierr)
    if(min_i1 /= expected_min ) then
       write(*,*) 'i1 min error'
       stop 1
    end if
    call adios2_variable_max(max_i1, vars_in(1), ierr)
    if(max_i1 /= expected_max ) then
       write(*,*) 'i1 max error'
       stop 1
    end if

    call adios2_inquire_variable(vars_in(2), ioGet, 'var_i2', ierr)
    if( vars_in(2)%valid .eqv. .false. ) then
       write(*,*) 'i2 invalid error'
       stop 1
    end if
    if( vars_in(2)%name /= 'var_i2' ) then
       write(*,*) 'i2 name error'
       stop 1
    end if
    if( vars_in(2)%type /= adios2_type_integer2 ) then
       write(*,*) 'i2 type error'
       stop 1
    end if
    if( vars_in(2)%ndims /= 3 ) then
       write(*,*) 'i2 dims error'
       stop 1
    end if
    call adios2_variable_min(min_i2, vars_in(2), ierr)
    if(min_i2 /= expected_min ) then
       write(*,*) 'i2 min error'
       stop 1
    end if
    call adios2_variable_max(max_i2, vars_in(2), ierr)
    if(max_i2 /= expected_max ) then
       write(*,*) 'i2 max error'
       stop 1
    end if

    call adios2_inquire_variable(vars_in(3), ioGet, 'var_i4', ierr)
    if( vars_in(3)%valid .eqv. .false. ) then
       write(*,*) 'i4 invalid error'
       stop 1
    end if
    if( vars_in(3)%name /= 'var_i4' ) then
       write(*,*) 'i4 name error'
       stop 1
    end if
    if( vars_in(3)%type /= adios2_type_integer4 ) then
       write(*,*) 'i4 type error'
       stop 1
    end if
    if( vars_in(3)%ndims /= 3 ) then
       write(*,*) 'i4 dims error'
       stop 1
    end if
    call adios2_variable_min(min_i4, vars_in(3), ierr)
    if(min_i4 /= expected_min ) then
       write(*,*) 'i4 min error'
       stop 1
    end if
    call adios2_variable_max(max_i4, vars_in(3), ierr)
    if(max_i4 /= expected_max ) then
       write(*,*) 'i4 max error'
       stop 1
    end if

    call adios2_inquire_variable(vars_in(4), ioGet, 'var_i8', ierr)
    if( vars_in(4)%valid .eqv. .false. ) then
       write(*,*) 'i8 invalid error'
       stop 1
    end if
    if( vars_in(4)%name /= 'var_i8' ) then
       write(*,*) 'i8 name error'
       stop 1
    end if
    if( vars_in(4)%type /= adios2_type_integer8 ) then
       write(*,*) 'i8 type error'
       stop 1
    end if
    if( vars_in(4)%ndims /= 3 ) then
       write(*,*) 'i8 dims error'
       stop 1
    end if
    call adios2_variable_min(min_i8, vars_in(4), ierr)
    if(min_i8 /= expected_min ) then
       write(*,*) 'i8 min error'
       stop 1
    end if
    call adios2_variable_max(max_i8, vars_in(4), ierr)
    if(max_i8 /= expected_max ) then
       write(*,*) 'i8 max error'
       stop 1
    end if

    call adios2_inquire_variable(vars_in(5), ioGet, 'var_r4', ierr)
    if( vars_in(5)%valid .eqv. .false. ) then
       write(*,*) 'r4 invalid error'
       stop 1
    end if
    if( vars_in(5)%name /= 'var_r4' ) then
       write(*,*) 'r4 name error'
       stop 1
    end if
    if( vars_in(5)%type /= adios2_type_real ) then
       write(*,*) 'r4 type error'
       stop 1
    end if
    if( vars_in(5)%ndims /= 3 ) then
       write(*,*) 'r4 dims error'
       stop 1
    end if
    call adios2_variable_min(min_r4, vars_in(5), ierr)
    if(min_r4 /= float(expected_min) ) then
       write(*,*) 'r4 min error'
       stop 1
    end if
    call adios2_variable_max(max_r4, vars_in(5), ierr)
    if(max_r4 /= float(expected_max) ) then
       write(*,*) 'r4 max error'
       stop 1
    end if

    call adios2_inquire_variable(vars_in(6), ioGet, 'var_r8', ierr)
    if( vars_in(6)%valid .eqv. .false. ) then
       write(*,*) 'r8 invalid error'
       stop 1
    end if
    if( vars_in(6)%name /= 'var_r8' ) then
       write(*,*) 'r8 name error'
       stop 1
    end if
    if( vars_in(6)%type /= adios2_type_dp ) then
       write(*,*) 'r8 type error'
       stop 1
    end if
    if( vars_in(6)%ndims /= 3 ) then
       write(*,*) 'r8 dims error'
       stop 1
    end if
    call adios2_variable_min(min_r8, vars_in(6), ierr)
    if(min_r8 /= float(expected_min) ) then
       write(*,*) 'r8 min error'
       stop 1
    end if
    call adios2_variable_max(max_r8, vars_in(6), ierr)
    if(max_r8 /= float(expected_max) ) then
       write(*,*) 'r8 max error'
       stop 1
    end if

    call adios2_get(bpReader, 'var_i1', in_data_i1, ierr)
    call adios2_get(bpReader, 'var_i2', in_data_i2, ierr)
    call adios2_get(bpReader, 'var_i4', in_data_i4, ierr)
    call adios2_get(bpReader, 'var_i8', in_data_i8, ierr)
    call adios2_get(bpReader, 'var_r4', in_data_r4, ierr)
    call adios2_get(bpReader, 'var_r8', in_data_r8, ierr)
    call adios2_end_step(bpReader, ierr)

    do k=1,isize*nz
      do j=1,ny
        do i=1,nx
           if(in_data_i1(i,j,k) /= current_step) then
              write(*,*) 'i1 read error'
              stop 1
           end if
           if(in_data_i2(i,j,k) /= current_step) then
              write(*,*) 'i2 read error'
              stop 1
           end if
           if(in_data_i4(i,j,k) /= current_step) then
              write(*,*) 'i4 read error'
              stop 1
           end if
           if(in_data_i8(i,j,k) /= current_step) then
              write(*,*) 'i8 read error'
              stop 1
           end if
           if(in_data_r4(i,j,k) /= REAL(current_step, 4)) then
              write(*,*) 'r4 read error'
              stop 1
           end if
           if(in_data_r8(i,j,k) /= current_step) then
              write(*,*) 'r8 read error'
              stop 1
           end if
        end do
      end do
    end do

  end do

  call adios2_close(bpReader, ierr)

  deallocate (in_data_i1)
  deallocate (in_data_i2)
  deallocate (in_data_i4)
  deallocate (in_data_i8)
  deallocate (in_data_r4)
  deallocate (in_data_r8)

  call adios2_finalize(adios, ierr)
  call MPI_Finalize(ierr)

end program TestBPWriteMemorySelectionRead3D
