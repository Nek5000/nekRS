 program TestBPWriteTypes
     use small_test_data
     use mpi
     use adios2
     implicit none

     integer(kind=8), dimension(1) :: count_dims, changing_count_dims
     integer :: inx, irank, isize, ierr, i, step_status, s

     type(adios2_adios) :: adios
     type(adios2_io) :: ioWrite, ioRead, ioDummy
     type(adios2_variable), dimension(13) :: variables
     type(adios2_engine) :: bpWriter, bpReader

     ! read handlers
     integer(kind=8) current_step, block_id

     ! data
     integer(kind=1), dimension(10) :: I8, inI8
     integer(kind=2), dimension(10) :: I16, inI16
     integer(kind=4), dimension(10) :: I32, inI32
     integer(kind=8), dimension(10) :: I64, inI64
     real(kind=4), dimension(10) :: R32, inR32
     real(kind=8), dimension(10) :: R64, inR64
     real(kind=8), dimension(10) :: changingR64, inchangingR64

     ! Program starts

     ! Launch MPI
     INTEGER provided

     ! MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
     call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
     call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

     ! Application variables
     inx = 10

     ! Variable dimensions
     count_dims(1) = inx
     changing_count_dims(1) = inx

     ! Create adios handler passing the communicator and error flag
     call adios2_init(adios, MPI_COMM_WORLD, ierr)

     !!!!!!!!!!!!!!!!!!!!!!!!! WRITER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Declare an IO process configuration inside adios
     call adios2_declare_io(ioWrite, adios, "ioWrite", ierr)

     ! Defines a variable to be written in bp format
     call adios2_define_variable(variables(1), ioWrite, "var_I8", &
                                 adios2_type_integer1, 1, &
                                 adios2_null_dims, adios2_null_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(2), ioWrite, "var_I16", &
                                 adios2_type_integer2, 1, &
                                 adios2_null_dims, adios2_null_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(3), ioWrite, "var_I32", &
                                 adios2_type_integer4, 1, &
                                 adios2_null_dims, adios2_null_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(4), ioWrite, "var_I64", &
                                 adios2_type_integer8, 1, &
                                 adios2_null_dims, adios2_null_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(5), ioWrite, "var_R32", &
                                 adios2_type_real, 1, &
                                 adios2_null_dims, adios2_null_dims, count_dims, &
                                 adios2_constant_dims, ierr)

     call adios2_define_variable(variables(6), ioWrite, "var_R64", &
                                 adios2_type_dp, 1, &
                                 adios2_null_dims, adios2_null_dims, count_dims, &
                                 .false., ierr)

     ! Global variables
     call adios2_define_variable(variables(7), ioWrite, "gvar_I8", &
                                 adios2_type_integer1, ierr)

     call adios2_define_variable(variables(8), ioWrite, "gvar_I16", &
                                 adios2_type_integer2, ierr)

     call adios2_define_variable(variables(9), ioWrite, "gvar_I32", &
                                 adios2_type_integer4, ierr)

     call adios2_define_variable(variables(10), ioWrite, "gvar_I64", &
                                 adios2_type_integer8, ierr)

     call adios2_define_variable(variables(11), ioWrite, "gvar_R32", &
                                 adios2_type_real, ierr)

     call adios2_define_variable(variables(12), ioWrite, "gvar_R64", &
                                 adios2_type_dp, ierr)

     call adios2_define_variable(variables(13), ioWrite, "var_changingR64", &
                                 adios2_type_dp, 1, &
                                 adios2_null_dims, adios2_null_dims, changing_count_dims, &
                                 adios2_variable_dims, ierr)

     write (*, *) "Engine type: ", ioWrite%engine_type
     if (TRIM(ioWrite%engine_type) /= 'File') then
        write(*,*) 'Wrong engine_type'
        stop 1
     end if

     call adios2_set_engine(ioWrite, "SST", ierr)
     write (*, *) "Engine type: ", ioWrite%engine_type
     if (TRIM(ioWrite%engine_type) /= 'SST') then
        write(*,*) 'Wrong engine_type'
        stop 1
     end if

     call adios2_at_io(ioDummy, adios, "ioWrite", ierr)
     write (*, *) "Engine type: ", ioDummy%engine_type
     if (TRIM(ioDummy%engine_type) /= 'SST') then
        write(*,*) 'Wrong engine_type'
        stop 1
     end if

     call adios2_set_engine(ioWrite, "BPFile", ierr)

     ! Open myVector_f.bp in write mode, this launches an engine
     call adios2_open(bpWriter, ioWrite, "ftypes_local.bp", adios2_mode_write, &
                      ierr)

     do s = 1, 3
         call adios2_begin_step(bpWriter, ierr)

         call adios2_current_step(current_step, bpWriter, ierr)
         if (current_step /= s - 1) then
            write(*,*) 'wrong current step'
            stop 1
         end if

         if (irank == 0 .and. s == 1) then
             call adios2_put(bpWriter, variables(7), data_I8(1), ierr)
             call adios2_put(bpWriter, variables(8), data_I16(1), ierr)
             call adios2_put(bpWriter, variables(9), data_I32(1), ierr)
             call adios2_put(bpWriter, variables(10), data_I64(1), ierr)
             call adios2_put(bpWriter, variables(11), data_R32(1), ierr)
             call adios2_put(bpWriter, variables(12), data_R64(1), ierr)
         end if

         do i = 1, inx
             I8(i) = data_I8(i) + INT(current_step, 1)
             I16(i) = data_I16(i) + INT(irank + current_step, 2)
             I32(i) = data_I32(i) + INT(irank + current_step, 4)
             I64(i) = data_I64(i) + INT(irank + current_step, 8)
             R32(i) = data_R32(i) + INT(irank + current_step, 4)
             R64(i) = data_R64(i) + INT(irank + current_step, 8)
         end do

         ! changing count
         changing_count_dims(1) = s
         call adios2_set_selection(variables(13), 1, adios2_null_dims, &
                                    changing_count_dims, ierr)
         do i = 1, s
             changingR64(i) = data_R64(i) + irank + current_step
         end do

         call adios2_put(bpWriter, variables(1), I8, ierr)
         call adios2_put(bpWriter, variables(2), I16, ierr)
         call adios2_put(bpWriter, variables(3), I32, ierr)
         call adios2_put(bpWriter, variables(4), I64, ierr)
         call adios2_put(bpWriter, variables(5), R32, ierr)
         call adios2_put(bpWriter, variables(6), R64, ierr)
         call adios2_put(bpWriter, variables(13), changingR64, ierr)

         call adios2_end_step(bpWriter, ierr)
     end do

     ! Closes engine1 and deallocates it, becomes unreachable
     call adios2_close(bpWriter, ierr)

     !!!!!!!!!!!!!!!!!!!!!!!!! READER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Declare io reader
     call adios2_declare_io(ioRead, adios, "ioRead", ierr)
     ! Open bpReader engine
     call adios2_open(bpReader, ioRead, "ftypes_local.bp", adios2_mode_read, ierr)

     do
         call adios2_begin_step(bpReader, adios2_step_mode_read, -1.0, &
                                step_status, ierr)
         if (step_status == adios2_step_status_end_of_stream) exit

         call adios2_current_step(current_step, bpReader, ierr)

         if (current_step == 0) then
             call adios2_inquire_variable(variables(7), ioRead, "gvar_I8", ierr)
             if (variables(7)%name /= 'gvar_I8') then
                write(*,*) 'gvar_I8 name not recognized'
                stop 1
             end if
             if (variables(7)%type /= adios2_type_integer1) then
                write(*,*) 'gvar_I8 type not recognized'
                stop 1
             end if

             call adios2_inquire_variable(variables(8), ioRead, "gvar_I16", ierr)
             if (variables(8)%name /= 'gvar_I16') then
                write(*,*) 'gvar_I16 name not recognized'
                stop 1
             end if
             if (variables(8)%type /= adios2_type_integer2) then
                write(*,*) 'gvar_I16 type not recognized'
                stop 1
             end if

             call adios2_inquire_variable(variables(9), ioRead, "gvar_I32", ierr)
             if (variables(9)%name /= 'gvar_I32') then
                write(*,*) 'gvar_I32 name not recognized'
                stop 1
             end if
             if (variables(9)%type /= adios2_type_integer4) then
                write(*,*) 'gvar_I32 type not recognized'
                stop 1
             end if

             call adios2_inquire_variable(variables(10), ioRead, "gvar_I64", ierr)
             if (variables(10)%name /= 'gvar_I64') then
                write(*,*) 'gvar_I64 name not recognized'
                stop 1
             end if
             if (variables(10)%type /= adios2_type_integer8) then
                write(*,*) 'gvar_I64 type not recognized'
                stop 1
             end if

             call adios2_inquire_variable(variables(11), ioRead, "gvar_R32", ierr)
             if (variables(11)%name /= 'gvar_R32') then
                write(*,*) 'gvar_R32 name not recognized'
                stop 1
             end if
             if (variables(11)%type /= adios2_type_real) then
                write(*,*) 'gvar_I64 type not recognized'
                stop 1
             end if

             call adios2_inquire_variable(variables(12), ioRead, "gvar_R64", ierr)
             if (variables(12)%name /= 'gvar_R64') then
                write(*,*) 'gvar_R64 name not recognized'
                stop 1
             end if
             if (variables(12)%type /= adios2_type_dp) then
                write(*,*) 'gvar_I64 type not recognized'
                stop 1
             end if
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

         call adios2_inquire_variable(variables(2), ioRead, "var_I16", ierr)
         if (variables(2)%name /= 'var_I16') then
            write(*,*) 'var_I16 not recognized'
            stop 1
         end if
         if (variables(2)%type /= adios2_type_integer2) then
            write(*,*) 'var_I16 type not recognized'
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

         call adios2_inquire_variable(variables(4), ioRead, "var_I64", ierr)
         if (variables(4)%name /= 'var_I64') then
            write(*,*) 'var_I64 not recognized'
            stop 1
         end if
         if (variables(4)%type /= adios2_type_integer8) then
            write(*,*) 'var_I64 type not recognized'
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

         call adios2_inquire_variable(variables(6), ioRead, "var_R64", ierr)
         if (variables(6)%name /= 'var_R64') then
            write(*,*) 'var_R64 not recognized'
            stop 1
         end if
         if (variables(6)%type /= adios2_type_dp) then
            write(*,*) 'var_R64 type not recognized'
            stop 1
         end if

         call adios2_inquire_variable(variables(13), ioRead, "var_changingR64", ierr)
         if (variables(13)%name /= 'var_changingR64') then
            write(*,*) 'var_changingR64 not recognized'
            stop 1
         end if
         if (variables(13)%type /= adios2_type_dp) then
            write(*,*) 'var_R64 type not recognized'
            stop 1
         end if

         do block_id = 0, isize - 1

             call adios2_set_block_selection(variables(1), block_id, ierr)
             call adios2_set_block_selection(variables(2), block_id, ierr)
             call adios2_set_block_selection(variables(3), block_id, ierr)
             call adios2_set_block_selection(variables(4), block_id, ierr)
             call adios2_set_block_selection(variables(5), block_id, ierr)
             call adios2_set_block_selection(variables(6), block_id, ierr)
             call adios2_set_block_selection(variables(13), block_id, ierr)

             call adios2_get(bpReader, variables(1), inI8, ierr)
             call adios2_get(bpReader, variables(2), inI16, ierr)
             call adios2_get(bpReader, variables(3), inI32, ierr)
             call adios2_get(bpReader, variables(4), inI64, ierr)
             call adios2_get(bpReader, variables(5), inR32, ierr)
             call adios2_get(bpReader, variables(6), inR64, ierr)
             call adios2_get(bpReader, variables(13), inchangingR64, ierr)
             call adios2_perform_gets(bpReader, ierr)

             do i = 1, inx
                 I8(i) = data_I8(i) + INT(current_step, 1)
                 I16(i) = data_I16(i) + INT(block_id + current_step, 2)
                 I32(i) = data_I32(i) + INT(block_id + current_step, 4)
                 I64(i) = data_I64(i) + INT(block_id + current_step, 8)
                 R32(i) = data_R32(i) + INT(block_id + current_step, 4)
                 R64(i) = data_R64(i) + INT(block_id + current_step, 8)
             end do

             do i = 1, inx
                if (I8(i) /= inI8(i)) then
                   write(*,*) 'Error reading var_I8'
                   stop 1
                end if
                if (I16(i) /= inI16(i)) then
                   write(*,*) 'Error reading var_I16'
                   stop 1
                end if
                if (I32(i) /= inI32(i)) then
                   write(*,*) 'Error reading var_I32'
                   stop 1
                end if
                if (I64(i) /= inI64(i)) then
                   write(*,*) 'Error reading var_I64'
                   stop 1
                end if
                if (R32(i) /= inR32(i)) then
                   write(*,*) 'Error reading var_R32'
                   stop 1
                end if
                if (R64(i) /= inR64(i)) then
                   write(*,*) 'Error reading var_R64'
                   stop 1
                end if

                 if( i < current_step) then
                    if (R64(i) /= inchangingR64(i)) then
                       write(*,*) 'Error reading var_changingR64'
                       stop 1
                    end if
                 end if

             end do

         end do

         call adios2_end_step(bpReader, ierr)
     end do

     call adios2_close(bpReader, ierr)

     ! Deallocates adios and calls its destructor
     call adios2_finalize(adios, ierr)

     call MPI_Finalize(ierr)

 end program TestBPWriteTypes
