from mpi4py import MPI
import numpy as np
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


sstReader = sstIO.Open("helloSst", adios2.Mode.Read)
loopStep = 0
while True:
    status = sstReader.BeginStep()
    if not status == adios2.StepStatus.OK:
        break
    var = sstIO.InquireVariable("bpFloats")
    shape = var.Shape()
    count = int(shape[0] / size)
    start = count * rank
    if rank == size - 1:
        count += shape[0] % size
    var.SetSelection([[start], [count]])
    # note: we need to use np.float32 to be compatible with data from C++ writer
    # using "float" works in Python only but leads to type mismatch with C++
    floatArray = np.zeros(count, dtype=np.float32)
    currentStep = sstReader.CurrentStep()
    sstReader.Get(var, floatArray)
    sstReader.EndStep()

    print(
        "Rank=",
        rank,
        "loop index =",
        loopStep,
        "stream step =",
        currentStep,
        "data =",
        floatArray,
        flush=True,
    )
    loopStep = loopStep + 1

print(f"Exited loop with StepStatus = {status}")
sstReader.Close()
