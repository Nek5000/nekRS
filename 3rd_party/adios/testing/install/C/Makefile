#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

ifeq ($(ADIOS2_HAVE_MPI),1)
	main_c=main_mpi.c
	mpiexec=$(MPIEXEC)
	test_mpi=test_mpi
else
	main_c=main_nompi.c
	mpiexec=
	test_mpi=
endif

all: test test_serial $(test_mpi)

test: adios_c_test adios_c_test_2
	$(mpiexec) ./adios_c_test
	$(mpiexec) ./adios_c_test_2

adios_c_test: main.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o adios_c_test main.o $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --c-libs)

main.o: $(main_c)
	$(CC) $(CFLAGS) $(ISYSROOT) -o main.o -c $(main_c) $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --c-flags)

adios_c_test_2: $(main_c)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(ISYSROOT) -o adios_c_test_2 $(main_c) $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -c)


test_serial: adios_c_serial_test adios_c_serial_test_2
	./adios_c_serial_test
	./adios_c_serial_test_2

adios_c_serial_test: main_nompi.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o adios_c_serial_test main_nompi.o $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --serial --c-libs)

main_nompi.o: main_nompi.c
	$(CC) $(CFLAGS) $(ISYSROOT) -o main_nompi.o -c main_nompi.c $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -s --c-flags)

adios_c_serial_test_2: main_nompi.c
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(ISYSROOT) -o adios_c_serial_test_2 main_nompi.c $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -s -c)


test_mpi: adios_c_mpi_test adios_c_mpi_test_2
	$(MPIEXEC) ./adios_c_mpi_test
	$(MPIEXEC) ./adios_c_mpi_test_2

adios_c_mpi_test: main_mpi.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o adios_c_mpi_test main_mpi.o $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --mpi --c-libs)

main_mpi.o: main_mpi.c
	$(CC) $(CFLAGS) $(ISYSROOT) -o main_mpi.o -c main_mpi.c $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -m --c-flags)

adios_c_mpi_test_2: main_mpi.c
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(ISYSROOT) -o adios_c_mpi_test_2 main_mpi.c $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -m -c)
