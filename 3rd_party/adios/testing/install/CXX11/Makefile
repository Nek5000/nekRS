#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

ifeq ($(ADIOS2_HAVE_MPI),1)
	main_cxx=main_mpi.cxx
	mpiexec=$(MPIEXEC)
	test_mpi=test_mpi
else
	main_cxx=main_nompi.cxx
	mpiexec=
	test_mpi=
endif

all: test test_serial $(test_mpi)

test: adios_cxx11_test adios_cxx11_test_2
	$(mpiexec) ./adios_cxx11_test
	$(mpiexec) ./adios_cxx11_test_2

adios_cxx11_test: main.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o adios_cxx11_test main.o $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --cxx-libs)

main.o: $(main_cxx)
	$(CXX) $(CXXFLAGS) $(ISYSROOT) -o main.o -c $(main_cxx) $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --cxx-flags)

adios_cxx11_test_2: $(main_cxx)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(ISYSROOT) -o adios_cxx11_test_2 $(main_cxx) $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -x)


test_serial: adios_cxx11_serial_test adios_cxx11_serial_test_2
	./adios_cxx11_serial_test
	./adios_cxx11_serial_test_2

adios_cxx11_serial_test: main_nompi.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o adios_cxx11_serial_test main_nompi.o $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --serial --cxx-libs)

main_nompi.o: main_nompi.cxx
	$(CXX) $(CXXFLAGS) $(ISYSROOT) -o main_nompi.o -c main_nompi.cxx $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -s --cxx-flags)

adios_cxx11_serial_test_2: main_nompi.cxx
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(ISYSROOT) -o adios_cxx11_serial_test_2 main_nompi.cxx $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -s -x)


test_mpi: adios_cxx11_mpi_test adios_cxx11_mpi_test_2
	$(MPIEXEC) ./adios_cxx11_mpi_test
	$(MPIEXEC) ./adios_cxx11_mpi_test_2

adios_cxx11_mpi_test: main_mpi.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o adios_cxx11_mpi_test main_mpi.o $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) --mpi --cxx-libs)

main_mpi.o: main_mpi.cxx
	$(CXX) $(CXXFLAGS) $(ISYSROOT) -o main_mpi.o -c main_mpi.cxx $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -m --cxx-flags)

adios_cxx11_mpi_test_2: main_mpi.cxx
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(ISYSROOT) -o adios_cxx11_mpi_test_2 main_mpi.cxx $(shell adios2-config$(ADIOS2_EXECUTABLE_SUFFIX) -m -x)
