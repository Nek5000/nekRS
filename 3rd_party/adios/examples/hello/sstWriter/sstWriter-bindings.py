from mpi4py import MPI
import numpy as np
from time import sleep
import adios2.bindings as adios2

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ADIOS MPI Communicator
adios = adios2.ADIOS(comm)

# ADIOS IO
sstIO = adios.DeclareIO("myIO")
sstIO.SetEngine("Sst")

# Python MPI is usually not compatible with MPI used to build ADIOS
# Must use TCP based data transport, or RDMA (ucx, libfabric) if available
sstIO.SetParameter("DataTransport", "WAN")

# note: we need to use np.float32 to be compatible with data from C++ writer
# using "float" works in Python only but leads to type mismatch with C++
myArray = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], dtype=np.float32)
myArray = 10.0 * rank + myArray
nx = len(myArray)
increment = nx * size * 1.0

# ADIOS Variable
ioArray = sstIO.DefineVariable(
    "bpFloats", myArray, [size * nx], [rank * nx], [nx], adios2.ConstantDims
)

sstWriter = sstIO.Open("helloSst", adios2.Mode.Write)
for i in range(4):
    print("Rank=", rank, "loop index =", i, "data =", myArray, flush=True)
    sstWriter.BeginStep()
    sstWriter.Put(ioArray, myArray, adios2.Mode.Sync)
    myArray += increment
    # Warning: the data is not published until EndStep is called
    sstWriter.EndStep()
    sleep(1.0)

sstWriter.Close()
