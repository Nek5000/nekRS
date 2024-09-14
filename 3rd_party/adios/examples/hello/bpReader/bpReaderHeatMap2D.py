#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# bpReaderHeatMap2D.py
#
#
#  Created on: Dec 5th, 2017
#      Author: William F Godoy godoywf@ornl.gov
#              Norbert Podhorszki pnorbert@ornl.gov
#

from mpi4py import MPI
import numpy
from adios2 import Stream, FileReader

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# User data
Nx = 10
Ny = 10

count = [Nx, Ny]
start = [rank * Nx, 0]
shape = [size * Nx, Ny]

temperatures = numpy.zeros(Nx * Ny, dtype=int)

for i in range(0, Nx):
    iGlobal = start[0] + i

    for j in range(0, Ny):
        value = iGlobal * shape[1] + j
        temperatures[i * Nx + j] = value

with Stream("HeatMap2D_py.bp", "w", comm) as obpStream:
    obpStream.write("temperature2D", temperatures, shape, start, count)
    if not rank:
        obpStream.write("N", [size, Nx, Ny])  # will be an array in output
        obpStream.write("Nx", numpy.array(Nx))  # will be a scalar in output
        obpStream.write("Ny", Ny)  # will be a scalar in output
        obpStream.write_attribute("nproc", size)  # will be a single value attribute in output
        obpStream.write_attribute("dimensions", [size * Nx, Ny], "temperature2D")

if not rank:
    with FileReader("HeatMap2D_py.bp", MPI.COMM_SELF) as ibpFile:
        # scalar variables are read as a numpy array with 0 dimension
        in_nx = ibpFile.read("Nx")
        in_ny = ibpFile.read("Ny")
        print(f"Incoming nx, ny = {in_nx}, {in_ny}")

        # single value attribute is read as a numpy array with 0 dimension
        in_nproc = ibpFile.read_attribute("nproc")
        print(f"Incoming nproc = {in_nproc}")
        # array attribute is read as a numpy array or string list
        in_dims = ibpFile.read_attribute("temperature2D/dimensions")
        print(f"Incoming diumensions = {in_dims}")

        # On option is to inquire a variable to know its type, shape
        # directly, not as strings, and then we can use the variable
        # object to set selection and/or set steps to read
        var_inTemperature = ibpFile.inquire_variable("temperature2D")
        if var_inTemperature is not None:
            var_inTemperature.set_selection([[2, 2], [4, 4]])
            inTemperatures = ibpFile.read(var_inTemperature)
            print(
                f"Incoming temperature map with selection "
                f"start = {var_inTemperature.start()}, count = {var_inTemperature.count()}"
            )
            for i in range(0, inTemperatures.shape[1]):
                print(str(inTemperatures[i]))
