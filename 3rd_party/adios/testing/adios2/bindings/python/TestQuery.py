#
from mpi4py import MPI
import numpy as np
import adios2.bindings as adios2
import sys

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# #######################################
# #  usage: <exe> [bp4 | bp5=default]  ##
# #######################################
numSteps = 5
queryFile = "query.xml"
targetVarName = "var0"

# User data
myArray = np.array([0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
Nx = myArray.size

# ADIOS MPI Communicator
adios = adios2.ADIOS(comm)

supportedEngines = ["bp5", "bp4"]
engineType = "bp5"
if len(sys.argv) > 1:
    engineType = sys.argv[1].lower()

if engineType in supportedEngines:
    if rank == 0:
        print("Using engine type:", engineType)
else:
    sys.exit("specified engine does not exist")


dataFileName = "test_" + engineType + ".bp"


def writeDataFile():
    bpIO = adios.DeclareIO("Writer")
    bpIO.SetEngine(engineType)

    ioArray = bpIO.DefineVariable(
        targetVarName, myArray, [size * Nx], [rank * Nx], [Nx], adios2.ConstantDims
    )

    bpFileWriter = bpIO.Open(dataFileName, adios2.Mode.Write)

    for i in range(numSteps):
        bpFileWriter.BeginStep()
        bpFileWriter.Put(ioArray, i * 10.0 + myArray / (rank + 1), adios2.Mode.Sync)
        bpFileWriter.EndStep()

    bpFileWriter.Close()


def createQueryFile():
    print(".. Writing query file to: ", queryFile)

    file1 = open(queryFile, "w")
    queryContent = [
        '<?xml version="1.0"?>\n',
        "<adios-query>\n",
        '  <io name="query">\n' '  <var name="' + targetVarName + '">\n',
        '    <op value="AND">\n',
        '      <range  compare="LT" value="15.0"/>\n',
        '      <range  compare="GT" value="4.0"/>\n',
        "    </op>\n",
        "  </var>\n",
        "  </io>\n",
        "</adios-query>\n",
    ]
    file1.writelines(queryContent)
    file1.close()


def doAnalysis(reader, touched_blocks, varList):
    print(" Step: ", reader.CurrentStep(), "  num touched blocks: ", len(touched_blocks))
    if 0 == reader.CurrentStep():
        assert len(touched_blocks) == min(size, 2)
    if 1 == reader.CurrentStep():
        assert len(touched_blocks) == size
    if 1 < reader.CurrentStep():
        assert len(touched_blocks) == 0

    values = []
    data = {}

    for var in varList:
        data[var] = []

    if len(touched_blocks) > 0:
        for n in touched_blocks:
            for var in varList:
                values = np.zeros(n[1], dtype=np.double)
                var.SetSelection(n)
                reader.Get(var, values, adios2.Mode.Sync)
                data[var].extend(values)


def queryDataFile():
    # # use no mpi
    adios_nompi = adios2.ADIOS()
    queryIO = adios_nompi.DeclareIO("query")

    reader = queryIO.Open(dataFileName, adios2.Mode.Read)
    print("dataFile=", dataFileName, "queryFile=", queryFile)
    touched_blocks = []

    print("Num steps: ", reader.Steps())

    while reader.BeginStep() == adios2.StepStatus.OK:
        # bp5 loads metadata after beginstep(),
        # therefore query has to be called per step
        w = adios2.Query(queryFile, reader)
        # assume only rank 0 wants to process result
        var = [queryIO.InquireVariable(targetVarName)]

        if rank == 0:
            touched_blocks = w.GetResult()
            doAnalysis(reader, touched_blocks, var)

        reader.EndStep()
    reader.Close()


def cleanUp():
    import os
    import shutil

    os.remove(queryFile)
    shutil.rmtree(dataFileName)
    print("  Cleanup generated files: ", queryFile, dataFileName)


#
# actual setup:
#


writeDataFile()

if 0 == rank:
    createQueryFile()
    queryDataFile()
    cleanUp()
