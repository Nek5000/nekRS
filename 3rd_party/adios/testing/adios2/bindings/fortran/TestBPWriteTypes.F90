program TestBPWriteTypes
   use small_test_data
#if ADIOS2_USE_MPI
   use mpi
#endif
   use adios2
   implicit none

   integer(kind=8), dimension(1) :: shape_dims, start_dims, count_dims
   integer(kind=4) :: inx, irank, isize, ierr, i, step_status
   integer(kind=8) :: nsteps

   type(adios2_adios) :: adios
   type(adios2_io) :: ioWrite, ioRead
   type(adios2_variable), dimension(14) :: variables
   type(adios2_derived_variable) :: derived_variable
   type(adios2_engine) :: bpWriter, bpReader
   character(len=15) :: inString
   character(len=:), allocatable :: varName, param_value
   character(len=:), allocatable :: engineType
   logical :: result

   ! read local value as global array
   integer(kind=4), dimension(:), allocatable :: inRanks

   ! read handlers
   integer(kind=4) :: ndims
   integer(kind=8), dimension(:), allocatable :: shape_in

   character(len=4096), dimension(:), allocatable :: varnamelist
   type(adios2_namestruct) :: namestruct

#if ADIOS2_USE_MPI
   ! Launch MPI
   INTEGER provided

   ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
   call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)
#else
   irank = 0
   isize = 1
#endif

   ! Application variables
   inx = 10

   ! Variable dimensions
   shape_dims(1) = isize*inx
   start_dims(1) = irank*inx
   count_dims(1) = inx

   if( adios%valid .eqv. .true. ) then
      write(*,*) 'Invalid adios default'
      stop 1
   end if
   if( ioWrite%valid .eqv. .true. ) then
      write(*,*) 'Invalid io default'
      stop 1
   end if

   do i=1,12
      if( variables(i)%valid .eqv. .true. ) then
         write(*,*) 'Invalid variables default'
         stop 1
      end if
   end do

   if( bpWriter%valid .eqv. .true. ) then
      write(*,*) 'Invalid engine default'
      stop 1
   end if


   ! Create adios handler passing the communicator and error flag
#if ADIOS2_USE_MPI
   call adios2_init(adios, MPI_COMM_WORLD, ierr)
#else
   call adios2_init(adios, ierr)
#endif
   if( adios%valid .eqv. .false. ) then
      write(*,*) 'Invalid adios2_init'
      stop 1
   end if

   !!!!!!!!!!!!!!!!!!!!!!!! WRITER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Declare an IO process configuration inside adios
   call adios2_declare_io(ioWrite, adios, "ioWrite", ierr)
   if( ioWrite%valid .eqv. .false. ) then
      write(*,*) 'Invalid adios2_declare_io'
      stop 1
   end if

   call adios2_at_io(ioWrite, adios, "ioWrite", ierr)
   if( ioWrite%valid .eqv. .false. ) then
      write(*,*) 'Invalid adios2_at_io'
      stop 1
   end if

   call adios2_in_config_file(result, ioWrite, ierr)
   if( result .eqv. .true. ) then
      write(*,*) 'Invalid ioWrite adios2_in_config_file'
      stop 1
   end if

   call adios2_set_engine(ioWrite, 'File', ierr)

   call adios2_set_parameter(ioWrite, 'ProfileUnits', 'Microseconds', ierr)

   call adios2_get_parameter(param_value, ioWrite, 'ProfileUnits', ierr)
   if( param_value /= "Microseconds") then
      write(*,*) 'Failed adios2_get_parameter ProfileUnits'
      stop 1
   end if

   call adios2_set_parameters(ioWrite, 'Threads=2, CollectiveMetadata = OFF', ierr)

   call adios2_get_parameter(param_value, ioWrite, 'Threads', ierr)
   if( param_value /= "2") then
      write(*,*) 'Failed adios2_get_parameter Threads'
      stop 1
   end if

   call adios2_get_parameter(param_value, ioWrite, 'CollectiveMetadata', ierr)
   if( param_value /= "OFF") then
      write(*,*) 'Failed adios2_get_parameter CollectiveMetadata'
      stop 1
   end if

   ! set back the default to make sure writing/reading test works
   call adios2_clear_parameters(ioWrite, ierr)
   call adios2_get_parameter(param_value, ioWrite, 'CollectiveMetadata', ierr)
   if( param_value /= "") then
      write(*,*) 'Still Could retrieve parameter CollectiveMetadata after clearing all parameters'
      stop 1
   end if

   deallocate(param_value)


   ! Defines a variable to be written in bp format
   call adios2_define_variable(variables(1), ioWrite, "var_I8", &
      adios2_type_integer1, 1, &
      shape_dims, start_dims, count_dims, &
      adios2_constant_dims, ierr)

   call adios2_define_variable(variables(2), ioWrite, "var_I16", &
      adios2_type_integer2, 1, &
      shape_dims, start_dims, count_dims, &
      adios2_constant_dims, ierr)

   call adios2_define_variable(variables(3), ioWrite, "var_I32", &
      adios2_type_integer4, 1, &
      shape_dims, start_dims, count_dims, &
      adios2_constant_dims, ierr)

   call adios2_define_variable(variables(4), ioWrite, "var_I64", &
      adios2_type_integer8, 1, &
      shape_dims, start_dims, count_dims, &
      adios2_constant_dims, ierr)

   call adios2_define_variable(variables(5), ioWrite, "var_R32", &
      adios2_type_real, 1, &
      shape_dims, start_dims, count_dims, &
      adios2_constant_dims, ierr)

   call adios2_define_variable(variables(6), ioWrite, "var_R64", &
      adios2_type_dp, 1, &
      shape_dims, start_dims, count_dims, &
      adios2_constant_dims, ierr)

   ! Global variables
   call adios2_define_variable(variables(7), ioWrite, "gvar_I8", &
      adios2_type_integer1,  ierr)

   call adios2_define_variable(variables(8), ioWrite, "gvar_I16", &
      adios2_type_integer2,  ierr)

   call adios2_define_variable(variables(9), ioWrite, "gvar_I32", &
      adios2_type_integer4,  ierr)

   call adios2_define_variable(variables(10), ioWrite, "gvar_I64", &
      adios2_type_integer8,  ierr)

   call adios2_define_variable(variables(11), ioWrite, "gvar_R32", &
      adios2_type_real,  ierr)

   call adios2_define_variable(variables(12), ioWrite, "gvar_R64", &
      adios2_type_dp,  ierr)

   call adios2_define_variable(variables(13), ioWrite, "gvar_Str", &
      adios2_type_string, ierr)

   ! local value
   call adios2_define_variable(variables(14), ioWrite, "lvar_i32", &
      adios2_type_integer4, &
      1, (/ adios2_local_value_dim /), &
      adios2_null_dims, &
      adios2_null_dims, &
      adios2_constant_dims, ierr)

   ! derived variable
   call adios2_define_derived_variable(derived_variable, ioWrite, "derived/magnitude_of_var_R64", &
      "x=var_R64 y=var_R64 z=var_R64 magnitude(x,y,z)", adios2_derived_var_type_metadata_only, ierr)
#if ADIOS2_HAVE_Derived_Variable
#define TOTAL_VAR_COUNT 15
#else
#define TOTAL_VAR_COUNT 14
#endif
   do i=1,13
      if( variables(i)%valid .eqv. .false. ) then
         write(*,*) 'Invalid adios2_define_variable'
         stop 1
      end if
   end do

   ! Testing adios2_variable_name for just two cases
   call adios2_variable_name(varName, variables(1), ierr)
   if (varName /= 'var_I8') then
      write(*,*) 'Invalid adios2_variable_name'
      stop 1
   end if

   call adios2_variable_name(varName, variables(2), ierr)
   if (varName /= 'var_I16') then
      write(*,*) 'Invalid adios2_variable_name'
      stop 1
   end if

   deallocate(varName)

   ! Open myVector_f.bp in write mode, this launches an engine
   if( ioWrite%valid .eqv. .false. ) then
      write(*,*) 'Invalid adios2_io'
      stop 1
   end if
   if( bpWriter%valid .eqv. .true. ) then
      write(*,*) 'Invalid adios2_engine pre-open'
      stop 1
   end if

#if ADIOS2_USE_MPI
   call adios2_open(bpWriter, ioWrite, "ftypes_mpi.bp", adios2_mode_write, ierr)
#else
   call adios2_open(bpWriter, ioWrite, "ftypes.bp", adios2_mode_write, ierr)
#endif
   if( bpWriter%valid .eqv. .false. ) then
      write(*,*) 'Invalid adios2_engine post-open'
      stop 1
   end if
#if ADIOS2_USE_MPI
   if( TRIM(bpWriter%name) /= "ftypes_mpi.bp") then
#else
   if( TRIM(bpWriter%name) /= "ftypes.bp") then
#endif
      write(*,*) 'Invalid adios2_engine name'
      stop 1
   end if

   if( TRIM(bpWriter%type) /= 'BP5Writer') then
      write(*,*) 'Engine Type ', TRIM(bpWriter%type)
      write(*,*) 'Invalid adios2_engine type'
      stop 1
   end if
   call adios2_io_engine_type(engineType, ioWrite, ierr)
   if( engineType /= 'File') then ! FIXME, different from the above!
      write(*,*) 'Engine Type ', engineType
      write(*,*) 'Invalid type from adios2_engine_type'
      stop 1
   end if

   if( bpWriter%mode /= adios2_mode_write) then
      write(*,*) 'Invalid adios2_engine mode'
      stop 1
   end if

   ! Put array contents to bp buffer, based on var1 metadata
   do i = 1, 3
      call adios2_begin_step(bpWriter, adios2_step_mode_append, -1.0, &
         step_status, ierr)

      if (irank == 0 .and. i == 1) then
         call adios2_put(bpWriter, variables(7), data_I8(1), ierr)
         call adios2_put(bpWriter, variables(8), data_I16(1), ierr)
         call adios2_put(bpWriter, variables(9), data_I32(1), ierr)
         call adios2_put(bpWriter, variables(10), data_I64(1), ierr)
         call adios2_put(bpWriter, variables(11), data_R32(1), ierr)
         call adios2_put(bpWriter, variables(12), data_R64(1), ierr)
         call adios2_put(bpWriter, variables(13), data_Strings(1), ierr)
      end if

      call adios2_put(bpWriter, variables(1), data_I8, ierr)
      call adios2_put(bpWriter, variables(2), data_I16, ierr)
      call adios2_put(bpWriter, variables(3), data_I32, ierr)
      call adios2_put(bpWriter, variables(4), data_I64, ierr)
      call adios2_put(bpWriter, variables(5), data_R32, ierr)
      call adios2_put(bpWriter, variables(6), data_R64, ierr)

      call adios2_put(bpWriter, variables(14), irank, ierr)

      call adios2_end_step(bpWriter, ierr)
   end do

   ! Closes engine1 and deallocates it, becomes unreachable
   call adios2_close(bpWriter, ierr)

   if( bpWriter%valid .eqv. .true. ) then
      write(*,*) 'Invalid adios2_close'
      stop 1
   end if

   !!!!!!!!!!!!!!!!!!!!!!!! READER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Declare io reader
   call adios2_declare_io(ioRead, adios, "ioRead", ierr)
   ! Open bpReader engine
#if ADIOS2_USE_MPI
   call adios2_open(bpReader, ioRead, "ftypes_mpi.bp", adios2_mode_readRandomAccess, ierr)
#else
   call adios2_open(bpReader, ioRead, "ftypes.bp", adios2_mode_readRandomAccess, ierr)
#endif
   call adios2_steps(nsteps, bpReader, ierr)
   if(nsteps /= 3) then
      write(*,*) 'ftypes.bp must have 3 steps'
      stop 1
   end if

   call adios2_available_variables(ioRead, namestruct, ierr)
   if (ierr /= 0) then
      write(*,*) 'adios2_available_variables returned with error'
      stop 1
   end if
   if (.not.namestruct%valid) then
      write(*,*) 'adios2_available_variables returned invalid struct'
      stop 1
   end if
   write(*,*) 'Number of variables = ', namestruct%count
   write(*,*) 'Max name length = ', namestruct%max_name_len
   if (namestruct%count /= TOTAL_VAR_COUNT) then
      write(*,*) 'adios2_available_variables returned not the expected 14'
      stop 1
   end if

   allocate(varnamelist(namestruct%count))

   call adios2_retrieve_names(namestruct, varnamelist, ierr)
   if (ierr /= 0) then
      write(*,*) 'adios2_retrieve_names returned with error'
      stop 1
   end if
   do i=1,namestruct%count
      write(*,'("Var[",i2,"] = ",a12)') i, varnamelist(i)
   end do
   deallocate(varnamelist)

   if (namestruct%f2c /= 0_8) then
      write(*,*) 'namestruct f2c pointer is not null after adios2_retrieve_names()'
      stop 1
   end if
   if (namestruct%valid) then
      write(*,*) 'namestruct is not invalidated after adios2_retrieve_names()'
      stop 1
   end if


   call adios2_inquire_variable(variables(1), ioRead, "var_I8", ierr)
   if (variables(1)%name /= 'var_I8') then
      write(*,*) 'var_I8 not recognized'
      stop 1
   end if
   if (variables(1)%type /= adios2_type_integer1) then
      write(*,*) 'var_I8 type not recognized'
      stop 1
   end if
   call adios2_variable_shape(shape_in, ndims, variables(1), ierr)
   if (ndims /= 1) then
      write(*,*) 'var_I8 ndims is not 1'
      stop 1
   end if
   if (shape_in(1) /= isize*inx) then
      write(*,*) 'var_I8 shape_in read failed'
      stop 1
   end if

   call adios2_inquire_variable(variables(2), ioRead, "var_I16", ierr)
   if (variables(2)%name /= 'var_I16') then
      write(*,*) 'var_I16 not recognized'
      stop 1
   end if
   if (variables(2)%type /= adios2_type_integer2) then
      write(*,*) 'var_I16 type not recognized'
      stop 1
   end if
   call adios2_variable_shape(shape_in, ndims, variables(2), ierr)
   if (ndims /= 1) then
      write(*,*) 'var_I16 ndims is not 1'
      stop 1
   end if
   if (shape_in(1) /= isize*inx) then
      write(*,*) 'var_I16 shape_in read failed'
      stop 1
   end if

   call adios2_inquire_variable(variables(3), ioRead, "var_I32", ierr)
   if (variables(3)%name /= 'var_I32') then
      write(*,*) 'var_I32 not recognized'
      stop 1
   end if
   if (variables(3)%type /= adios2_type_integer4) then
      write(*,*) 'var_I32 type not recognized'
      stop 1
   end if
   call adios2_variable_shape(shape_in, ndims, variables(3), ierr)
   if (ndims /= 1) then
      write(*,*) 'var_I32 ndims is not 1'
      stop 1
   end if
   if (shape_in(1) /= isize*inx) then
      write(*,*) 'var_I32 shape_in read failed'
      stop 1
   end if

   call adios2_inquire_variable(variables(4), ioRead, "var_I64", ierr)
   if (variables(4)%name /= 'var_I64') then
      write(*,*) 'var_I64 not recognized'
      stop 1
   end if
   if (variables(4)%type /= adios2_type_integer8) then
      write(*,*) 'var_I64 type not recognized'
      stop 1
   end if
   call adios2_variable_shape(shape_in, ndims, variables(4), ierr)
   if (ndims /= 1) then
      write(*,*) 'var_I64 ndims is not 1'
      stop 1
   end if
   if (shape_in(1) /= isize*inx) then
      write(*,*) 'var_I64 shape_in read failed'
      stop 1
   end if

   call adios2_inquire_variable(variables(5), ioRead, "var_R32", ierr)
   if (variables(5)%name /= 'var_R32') then
      write(*,*) 'var_R32 not recognized'
      stop 1
   end if
   if (variables(5)%type /= adios2_type_real) then
      write(*,*) 'var_R32 type not recognized'
      stop 1
   end if
   call adios2_variable_shape(shape_in, ndims, variables(5), ierr)
   if (ndims /= 1) then
      write(*,*) 'var_R32 ndims is not 1'
      stop 1
   end if
   if (shape_in(1) /= isize*inx) then
      write(*,*) 'var_R32 shape_in read failed'
      stop 1
   end if

   call adios2_inquire_variable(variables(6), ioRead, "var_R64", ierr)
   if (variables(6)%name /= 'var_R64') then
      write(*,*) 'var_R64 not recognized'
      stop 1
   end if
   if (variables(6)%type /= adios2_type_dp) then
      write(*,*) 'var_R64 type not recognized'
      stop 1
   end if
   call adios2_variable_shape(shape_in, ndims, variables(6), ierr)
   if (ndims /= 1) then
      write(*,*) 'var_R64 ndims is not 1'
      stop 1
   end if
   if (shape_in(1) /= isize*inx) then
      write(*,*) 'var_R64 shape_in read failed'
      stop 1
   end if

   call adios2_inquire_variable(variables(13), ioRead, "gvar_Str", ierr)
   call adios2_get(bpReader, variables(13), inString, ierr)
   call adios2_perform_gets(bpReader, ierr)
   if( inString /= data_Strings(1) ) then
      write(*,*) 'gvar_Str read failed'
      stop 1
   end if

   call adios2_inquire_variable(variables(14), ioRead, "lvar_i32", ierr)
   allocate(inRanks(isize))
   call adios2_get(bpReader, variables(14), inRanks, ierr)
   call adios2_perform_gets(bpReader, ierr)
   if( inRanks(irank+1) /= irank ) then
      write(*,*) 'lvar_i32 read failed'
      stop 1
   end if
   deallocate(inRanks)

   call adios2_close(bpReader, ierr)

   ! Deallocates adios and calls its destructor
   call adios2_finalize(adios, ierr)
   if( adios%valid .eqv. .true. ) stop 'Invalid adios2_finalize'

#if ADIOS2_USE_MPI
   call MPI_Finalize(ierr)
#endif

end program TestBPWriteTypes
