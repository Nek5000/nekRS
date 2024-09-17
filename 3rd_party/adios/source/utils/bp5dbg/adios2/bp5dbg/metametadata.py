from os import fstat

import numpy as np


def ReadMetaMetadataRecord(buf, rec, pos, fileSize):
    # Read one metametadata record
    startoffset = str(pos).rjust(8)
    # 8 bytes MetaMetaIDLen + MetaMetaInfoLen Length
    mmIDlen = np.frombuffer(buf, dtype=np.uint64, count=1, offset=pos)
    pos = pos + 8
    mmInfolen = np.frombuffer(buf, dtype=np.uint64, count=1, offset=pos)
    pos = pos + 8
    # mmid = np.frombuffer(buf, dtype=np.uint8, count=mmIDlen[0], offset=pos)
    pos = pos + int(mmIDlen[0])
    # mminfo = np.frombuffer(buf, dtype=np.uint8,
    #                        count=mmInfolen[0], offset=pos)
    pos = pos + int(mmInfolen[0])

    recs = str(rec).rjust(7)
    idlen = str(mmIDlen[0]).rjust(10)
    infolen = str(mmInfolen[0]).rjust(12)

    print(
        f" | {recs} | {startoffset} | {idlen} | {infolen} |")

    return pos


def DumpMetaMetaData(fileName):
    print("========================================================")
    print("    MetaMetadata File: " + fileName)
    print("========================================================")
    with open(fileName, "rb") as f:
        fileSize = fstat(f.fileno()).st_size
        print(f" File size = {fileSize}")
        buf = f.read(fileSize)
        print(" --------------------------------------------------")
        print(" |  Record |  Offset  |  ID length |  Info length |")
        print(" --------------------------------------------------")
        pos = 0
        rec = 0
        while (pos < fileSize - 1):
            pos = ReadMetaMetadataRecord(buf, rec, pos, fileSize)
            rec = rec + 1
        print(" --------------------------------------------------")
    return True


if __name__ == "__main__":
    print("ERROR: Utility main program is bp5dbg.py")
