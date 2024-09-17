#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

ifeq ($(ADIOS2_HAVE_MPI),1)
	main_f90=main_mpi.f90
	mpiexec=$(MPIEXEC)
	test_mpi=test_mpi
else
	main_f90=main_nompi.f90
	mpiexec=
	test_mpi=
endif

all: test test_serial $(test_mpi)

test: adios_fortran_test adios_fortran_test_2
	$(mpiexec) ./adios_fortran_test
	$(mpiexec) ./adios_fortran_test_2

adios_fortran_test: main.o
	$(FC) $(FFLAGS) -o adios_fortran_test main.o $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --fortran-libs)

main.o: $(main_f90)
	$(FC) $(FFLAGS) -o main.o -c $(main_f90) $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --fortran-flags)

adios_fortran_test_2: $(main_f90)
	$(FC) $(FFLAGS) -o adios_fortran_test_2 $(main_f90) $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -f)


test_serial: adios_fortran_serial_test adios_fortran_serial_test_2
	./adios_fortran_serial_test
	./adios_fortran_serial_test_2

adios_fortran_serial_test: main_nompi.o
	$(FC) $(FFLAGS) -o adios_fortran_serial_test main_nompi.o $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --serial --fortran-libs)

main_nompi.o: main_nompi.f90
	$(FC) $(FFLAGS) -o main_nompi.o -c main_nompi.f90 $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -s --fortran-flags)

adios_fortran_serial_test_2: main_nompi.f90
	$(FC) $(FFLAGS) -o adios_fortran_serial_test_2 main_nompi.f90 $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -s -f)


test_mpi: adios_fortran_mpi_test adios_fortran_mpi_test_2
	$(MPIEXEC) ./adios_fortran_mpi_test
	$(MPIEXEC) ./adios_fortran_mpi_test_2

adios_fortran_mpi_test: main_mpi.o
	$(FC) $(FFLAGS) -o adios_fortran_mpi_test main_mpi.o $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --mpi --fortran-libs)

main_mpi.o: main_mpi.f90
	$(FC) $(FFLAGS) -o main_mpi.o -c main_mpi.f90 $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -m --fortran-flags)

adios_fortran_mpi_test_2: main_mpi.f90
	$(FC) $(FFLAGS) -o adios_fortran_mpi_test_2 main_mpi.f90 $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -m -f)
