from os import fstat

import numpy as np


def ReadMetadataStep(f, fileSize, MetadataEntry, WriterMapEntry):
    # Read metadata of one step
    step = MetadataEntry['step']
    mdpos = MetadataEntry['mdpos']
    mdsize = MetadataEntry['mdsize']
    # flushcount = MetadataEntry['flushcount']
    WriterCount = WriterMapEntry['WriterCount']

    if mdpos + mdsize > fileSize:
        print(f"ERROR: step {step} metadata pos {mdpos} + size {mdsize} "
              f"is beyond the metadata file size {fileSize}")
        return False

    currentpos = f.tell()
    if mdpos > currentpos:
        print(f"Offset {currentpos}..{mdpos-1} is a gap unaccounted for")

    if mdpos < currentpos:
        print(f"ERROR: step {step} metadata pos {mdpos} points before the "
              f"expected position in file {currentpos}")
        return False

    f.seek(mdpos)
    buf = f.read(mdsize)
    pos = 0

    if step > 0:
        print("========================================================")
    print(f"Step {step}: ")
    print(f"  Offset = {mdpos}")

    mdsize_in_file = np.frombuffer(buf, dtype=np.uint64, count=1,
                                   offset=pos)
    pos = pos + 8
    if (mdsize == mdsize_in_file[0] + 8):
        print(f"  Size = {mdsize_in_file[0]}")
    else:
        print(f"ERROR: md record supposed to be {mdsize-8} + 8 bytes "
              f"(as recorded in index), but found in file "
              f"{mdsize_in_file[0]}")

    MDPosition = mdpos + 2 * 8 * WriterCount
    print("  Variable metadata entries: ")
    for w in range(0, WriterCount):
        a = np.frombuffer(buf, dtype=np.uint64, count=1, offset=pos)
        thisMDSize = int(a[0])
        pos = pos + 8
        print(f"    Writer {w}: md size {thisMDSize} "
              f"offset {MDPosition}")
        MDPosition = MDPosition + thisMDSize

    print("  Attribute metadata entries: ")
    for w in range(0, WriterCount):
        a = np.frombuffer(buf, dtype=np.uint64, count=1, offset=pos)
        thisMDSize = int(a[0])
        pos = pos + 8
        print(f"    Writer {w}: md size {thisMDSize} "
              f"offset {MDPosition}")
        MDPosition = MDPosition + thisMDSize

    if (mdsize_in_file != MDPosition - mdpos):
        print(f"ERROR: entries supposed to end at start offset+size "
              f"{mdpos}+{mdsize_in_file[0]}, but it ends instead on offset "
              f"{MDPosition}")

    return True


def DumpMetaData(fileName, MetadataIndexTable, WriterMap):
    print("========================================================")
    print("    Metadata File: " + fileName)
    print("========================================================")

    # print(f"MetadataIndexTable={MetadataIndexTable}")
    # print(f"WriterMap={WriterMap}")

    with open(fileName, "rb") as f:
        fileSize = fstat(f.fileno()).st_size
        for MetadataEntry in MetadataIndexTable:
            WriterMapEntry = WriterMap[MetadataEntry['writermapindex']]
            status = ReadMetadataStep(
                f, fileSize, MetadataEntry, WriterMapEntry)
    return status


if __name__ == "__main__":
    print("ERROR: Utility main program is bp5dbg.py")
