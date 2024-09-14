#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# bpReaderHeatMap2D-bindings.py
#
#
#  Created on: Dec 5th, 2017
#      Author: William F Godoy godoywf@ornl.gov
#              Norbert Podhorszki pnorbert@ornl.gov
#

from mpi4py import MPI
import numpy
import adios2.bindings as adios2

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
#        print(str(i) + "," + str(j) + " " + str(value))


# ADIOS portion
adios = adios2.ADIOS(comm)
ioWrite = adios.DeclareIO("ioWriter")

varTemperature = ioWrite.DefineVariable(
    "temperature2D", temperatures, shape, start, count, adios2.ConstantDims
)

obpStream = ioWrite.Open("HeatMap2D_py_bindings.bp", adios2.Mode.Write)
obpStream.BeginStep()
obpStream.Put(varTemperature, temperatures)
obpStream.EndStep()
obpStream.Close()


if rank == 0:
    ioRead = adios.DeclareIO("ioReader")

    ibpStream = ioRead.Open("HeatMap2D_py_bindings.bp", adios2.Mode.Read, MPI.COMM_SELF)

    ibpStream.BeginStep()

    var_inTemperature = ioRead.InquireVariable("temperature2D")

    if var_inTemperature is not None:
        var_inTemperature.SetSelection([[2, 2], [4, 4]])

        inSize = var_inTemperature.SelectionSize()
        print("Incoming size " + str(inSize))
        inTemperatures = numpy.zeros(var_inTemperature.Count(), dtype=int)

        ibpStream.Get(var_inTemperature, inTemperatures, adios2.Mode.Sync)

        print("Incoming temperature map")
        for i in range(0, inTemperatures.shape[1]):
            print(str(inTemperatures[i]) + " ")

    ibpStream.EndStep()
    ibpStream.Close()
