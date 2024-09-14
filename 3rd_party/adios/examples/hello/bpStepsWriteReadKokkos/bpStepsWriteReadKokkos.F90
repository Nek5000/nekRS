program TestBPWriteReadHeatMap2D use mpi use adios2

    implicit none

        integer(kind = 8)::sum_i1,
    sum_i2 type(adios2_adios)::adios type(adios2_io)::ioPut, ioGet type(adios2_engine)::bpWriter,
    bpReader type(adios2_variable), dimension(1)::var_g,
    var_gIn

    integer(kind = 2),
    dimension(
        :,
        :),
    allocatable ::g, &sel_g integer(kind = 8),
    dimension(2)::ishape, istart, icount integer(kind = 8),
    dimension(2)::sel_start, sel_count integer ::ierr, irank, isize,
    step_status integer ::in1, in2 integer ::i1,
    i2

    call MPI_INIT(ierr) call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr) call
    MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)

        in1 = 3 in2 = 4 icount = (/ in1, in2 /) istart = (/ 0, in2 *irank /) ishape = (/ in1,
                                                                                       in2 *isize /)

            allocate(g(in1, in2)) do i2 = 1,
        in2 do i1 = 1,
        in1 g(i1, i2) = irank + i1 end do end do

                        !Start adios2 Writer call adios2_init(adios, MPI_COMM_WORLD, ierr) call
                        adios2_declare_io(ioPut, adios, 'WriteIO', ierr)

                            call adios2_define_variable(var_g(1), ioPut, &'bpFloats',
                                                        adios2_type_integer2, &2, ishape, istart,
                                                        icount, &adios2_constant_dims, ierr)

                                call
                        adios2_open(bpWriter, ioPut, 'BPFortranKokkos.bp', adios2_mode_write, &ierr)

                            call adios2_put(bpWriter, var_g(1), g, ierr)

                                call adios2_close(bpWriter, ierr)

                                    if (allocated(g)) deallocate(g)

                                        !Start adios2 Reader in rank 0 if (irank == 0) then

                        call adios2_declare_io(ioGet, adios, 'ReadIO', ierr)

                            call adios2_open(bpReader, ioGet, 'BPFortranKokkos.bp',
                                             &adios2_mode_read, MPI_COMM_SELF, ierr)

                                call
                        adios2_begin_step(bpReader, adios2_step_mode_read, -1., &step_status, ierr)

                            call adios2_inquire_variable(var_gIn(1), ioGet, &'bpFloats', ierr)

                                sel_start = (/ 0, 0 /) sel_count = (/ ishape(1), ishape(2) /)

                            allocate(sel_g(ishape(1), ishape(2))) sel_g = 0

        call adios2_set_selection(var_gIn(1), 2, sel_start, sel_count, &ierr) call
        adios2_get(bpReader, var_gIn(1), sel_g, ierr)

            call adios2_end_step(bpReader, ierr)

                call adios2_close(bpReader, ierr)

                    do i2 = 1,
        INT(sel_count(2), 4) do i1 = 1,
        INT(sel_count(1), 4) write(6, "(i8)", advance = "no") sel_g(i1, i2) end do write(6, *) end
    do

    if (allocated(sel_g)) deallocate(sel_g)

        end if

    call adios2_finalize(adios, ierr) call MPI_Finalize(ierr)

        end program TestBPWriteReadHeatMap2D
