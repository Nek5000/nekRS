#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# bpWriter-bindings.py : Low-level Python API example
#  Created on: Feb 2, 2017
#      Authors: William F Godoy godoywf@ornl.gov
#               Norbert Podhorszki pnorbert@ornl.gov
#
# adios2.bindings gives access to the low level classes that are compiled from the C++ classes
# Only necessary to use this when the adios2 python API is insufficient for some reason
# This example doesn't use anything that is missing from the python API though.

from mpi4py import MPI
import numpy
import adios2.bindings as adios2

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# User data
myArray = numpy.array([0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
myArray = myArray + (rank * 10)
Nx = myArray.size

# ADIOS MPI Communicator
adios = adios2.ADIOS(comm)

# ADIOS IO
bpIO = adios.DeclareIO("BPFile_N2N")
bpIO.SetEngine("BPFile")
# bpIO.SetParameters( {"Threads" : "2", "ProfileUnits" : "Microseconds",
# "InitialBufferSize" : "17Kb"} )
bpIOParams = {}
bpIOParams["Threads"] = "2"
bpIOParams["ProfileUnits"] = "Microseconds"
bpIOParams["InitialBufferSize"] = "17Kb"
bpIO.SetParameters(bpIOParams)

fileID = bpIO.AddTransport("File", {"Library": "fstream"})

# ADIOS Variable name, shape, start, offset, constant dims
ioArray = bpIO.DefineVariable(
    "bpArray", myArray, [size * Nx], [rank * Nx], [Nx], adios2.ConstantDims
)

varNx = bpIO.DefineVariable("Nx", numpy.array(Nx))  # type is derived from numpy array type
bpIO.DefineAttribute("size", Nx, "bpArray")
bpIO.DefineAttribute("dimensions", ["Nx"], "bpArray")

# ADIOS Engine
bpFileWriter = bpIO.Open("bpWriter-py-bindings.bp", adios2.Mode.Write)
bpFileWriter.BeginStep()
bpFileWriter.Put(ioArray, myArray, adios2.Mode.Sync)
bpFileWriter.Put(varNx, numpy.array(Nx), adios2.Mode.Sync)
bpFileWriter.EndStep()
bpFileWriter.Close()

# Read content:
# bpls -la bpWriter-py-bindings.bp
# bpls -la bpWriter-py-bindings.bp -d bpArray -n 10
# bpls -la bpWriter-py-bindings.bp -d bpArray -n 10 -D
