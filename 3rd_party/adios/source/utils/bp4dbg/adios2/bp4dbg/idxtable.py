import numpy as np
from os import fstat
from .utils import *


def ReadIndex(f, fileSize):
    nBytes = fileSize - f.tell()
    if nBytes <= 0:
        return True
    nRows = int(nBytes / 64)
    table = f.read(nBytes)
    print(" ")
    print("-----------------------------------------------------------------"
          "------------------------------------------")
    print("|   Step   |   Rank   |    PGPtr    |   VarPtr    |   AttPtr    |"
          "   EndPtr    |  Timestamp  |    unused   |")
    print("+----------------------------------------------------------------"
          "-----------------------------------------+")
    for r in range(0, nRows):
        pos = r * 64
        data = np.frombuffer(table, dtype=np.uint64, count=8, offset=pos)
        step = str(data[0]).rjust(9)
        rank = str(data[1]).rjust(9)
        pgptr = str(data[2]).rjust(12)
        varptr = str(data[3]).rjust(12)
        attptr = str(data[4]).rjust(12)
        endptr = str(data[5]).rjust(12)
        time = str(data[6]).rjust(12)
        unused = str(data[7]).rjust(12)
        print("|" + step + " |" + rank + " |" + pgptr + " |" + varptr + " |" +
              attptr + " |" + endptr + " |" + time + " |" + unused + " |")

    print("-----------------------------------------------------------------"
          "------------------------------------------")

    if fileSize - f.tell() > 1:
        print("ERROR: There are {0} bytes at the end of file"
              " that cannot be interpreted".format(fileSize - f.tell() - 1))
        return False

    return True


def DumpIndexTable(fileName):
    print("========================================================")
    print("    Index Table File: " + fileName)
    print("========================================================")
    status = False
    with open(fileName, "rb") as f:
        fileSize = fstat(f.fileno()).st_size
        status = ReadHeader(f, fileSize, "Index Table")
        if status:
            status = ReadIndex(f, fileSize)
    return status


if __name__ == "__main__":
    print("ERROR: Utility main program is bp4dbg.py")
