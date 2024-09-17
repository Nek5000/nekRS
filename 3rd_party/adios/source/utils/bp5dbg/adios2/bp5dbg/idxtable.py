from os import fstat

import numpy as np

from .utils import ReadHeader

# metadata index table (list of dictionary)
#   step:           step (superfluous, since list idx = step)
#   mdpos:          pos in metadata file for step
#   mdsize:         size of metadata
#   flush_count:    flush count
#   writermapindex: entry in WriterMap relevant for this step

# WriterMap table (list of dictionary)
#   StartStep:       this entry is valid from this step
#   WriterCount:     number of writers
#   AggregatorCount: number of aggregators
#   SubfilesCount:   number of subfiles


def ReadWriterMap(bytearray, pos, verbose):
    data = np.frombuffer(bytearray, dtype=np.uint64, count=3,
                         offset=pos)
    WriterCount = int(data[0])
    AggregatorCount = int(data[1])
    SubfileCount = int(data[2])
    pos = pos + 3 * 8

    data = np.frombuffer(bytearray, dtype=np.uint64, count=WriterCount,
                         offset=pos)
    if verbose:
        print("  WriterMap: Writers = {0}  Aggregators = {1}  Subfiles = {2}"
              .format(WriterCount, AggregatorCount, SubfileCount))
        print("  =====================")
        print("  |  Rank  |  Subfile |")
        print("  ---------------------")
        for r in range(0, WriterCount):
            rank = str(r).rjust(7)
            sub = str(data[r]).rjust(8)
            print("  |" + rank + " | " + sub + " |")
        print("  =====================")

    pos = pos + WriterCount * 8
    return pos, WriterCount, AggregatorCount, SubfileCount


def ReadIndex(f, fileSize, verbose):
    nBytes = fileSize - f.tell()
    if nBytes <= 0:
        return True
    table = f.read(nBytes)
    WriterCount = 0
    pos = 0
    step = 0
    MetadataIndexTable = []
    WriterMap = []
    WriterMapIdx = -1
    while pos < nBytes:
        if verbose:
            print("-----------------------------------------------" +
                  "---------------------------------------------------")
        record = chr(table[pos])
        pos = pos + 1
        reclen = np.frombuffer(table, dtype=np.uint64, count=1,
                               offset=pos)
        pos = pos + 8
        if verbose:
            print("Record '{0}', length = {1}".format(record, reclen))
        if record == 's':
            # print("Step record, length = {0}".format(reclen))
            data = np.frombuffer(table, dtype=np.uint64, count=3,
                                 offset=pos)
            stepstr = str(step).ljust(6)
            mdatapos = str(data[0]).ljust(10)
            mdatasize = str(data[1]).ljust(10)
            flushcount = str(data[2]).ljust(3)
            FlushCount = data[2]

            md = {"step": step,
                  "mdpos": int(data[0]),
                  "mdsize": int(data[1]),
                  "flushcount": int(data[2]),
                  "writermapindex": WriterMapIdx}
            MetadataIndexTable.append(md)

            if verbose:
                print("|   Step = " + stepstr + "| MetadataPos = " +
                      mdatapos + " |  MetadataSize = " + mdatasize +
                      "   | FlushCount = " + flushcount + "|")

            pos = pos + 3 * 8

            for Writer in range(0, WriterCount):
                start = " Writer " + str(Writer) + " data "
                thiswriter = np.frombuffer(table, dtype=np.uint64,
                                           count=int(FlushCount * 2 + 1),
                                           offset=pos)

                for i in range(0, FlushCount):  # for flushes
                    start += ("loc:" + str(thiswriter[int(i * 2)]) + " siz:" +
                              str(thiswriter[i * 2 + 1]) + "; ")
                start += "loc:" + str(thiswriter[int(FlushCount * 2)])
                pos = int(pos + (FlushCount * 2 + 1) * 8)
                if verbose:
                    print(start)

            step = step + 1

        elif record == 'w':
            # print("WriterMap record")
            pos, WriterCount, AggregatorCount, SubfileCount = ReadWriterMap(
                table, pos, verbose)
            wmd = {"StartStep": step,
                   "WriterCount": WriterCount,
                   "AggregatorCount": AggregatorCount,
                   "SubfilesCount": SubfileCount}
            WriterMap.append(wmd)
            WriterMapIdx = WriterMapIdx + 1

        elif record == 'm':
            if verbose:
                print("MetaMeta record")

        else:
            if verbose:
                print("Unknown record {0}, lenght = {1}".format(
                    record, reclen))

        if verbose:
            print("---------------------------------------------------" +
                  "-----------------------------------------------")

    if fileSize - f.tell() > 1:
        print("ERROR: There are {0} bytes at the end of file"
              " that cannot be interpreted".format(fileSize - f.tell() - 1))
        return False, MetadataIndexTable, WriterMap

    return True, MetadataIndexTable, WriterMap


def DumpIndexTable(fileName, verbose):
    if verbose:
        print("========================================================")
        print("    Index Table File: " + fileName)
        print("========================================================")
    status = False
    MetadataIndexTable = []
    WriterMap = []
    with open(fileName, "rb") as f:
        fileSize = fstat(f.fileno()).st_size
        status = ReadHeader(
            f, fileSize, "Index Table", verbose)
        if isinstance(status, list):
            status = status[0]
        if status:
            status, MetadataIndexTable, WriterMap = ReadIndex(
                f, fileSize, verbose)
    return status, MetadataIndexTable, WriterMap


if __name__ == "__main__":
    print("ERROR: Utility main program is bp5dbg.py")
