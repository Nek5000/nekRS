#!/usr/bin/env python3

import argparse
import glob
from os.path import basename, exists, isdir

from adios2.bp5dbg import DumpIndexTable, DumpMetaData, DumpMetaMetaData

MetadataIndexTable = []
WriterMap = []
status = True


def SetupArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "FILE", help="Name of the input file (.bp, .bp/md.idx, " +
        ".bp/md.0 or .bp/data.XXX)")
    # parser.add_argument("--printdata", "-p",
    #                    help="Dump data of this variable as well", default="")
    parser.add_argument("--verbose", "-v",
                        help="More verbosity", action="count")
    parser.add_argument("--no-indextable", "-x",
                        help="Do not print index table md.idx",
                        action="store_true")
    parser.add_argument("--no-metadata", "-m",
                        help="Do not print metadata md.0",
                        action="store_true")
    parser.add_argument("--no-metametadata", "-M",
                        help="Do not print meta-metadata mmd.0",
                        action="store_true")
    parser.add_argument("--no-data", "-d",
                        help="Do not print data data.*", action="store_true")
    args = parser.parse_args()

    # default values
    args.idxFileName = ""
    args.dumpIdx = False
    args.metametadataFileName = ""
    args.dumpMetaMetadata = False
    args.metadataFileName = ""
    args.dumpMetadata = False
    args.dataFileName = ""
    args.dumpData = False

    # print("Verbosity = {0}".format(args.verbose))
    return args


def CheckFileName(args):
    if not exists(args.FILE):
        print("ERROR: File " + args.FILE + " does not exist", flush=True)
        exit(1)
    if isdir(args.FILE):
        if not args.no_indextable:
            args.idxFileName = args.FILE + "/" + "md.idx"
            args.dumpIdx = True
        if not args.no_metadata:
            args.metadataFileName = args.FILE + "/" + "md.[0-9]*"
            args.dumpMetadata = True
        if not args.no_metametadata:
            args.metametadataFileName = args.FILE + "/" + "mmd.[0-9]*"
            args.dumpMetaMetadata = True
        if not args.no_data:
            args.dataFileName = args.FILE + "/" + "data.[0-9]*"
            args.dumpData = True
        return

    name = basename(args.FILE)
    if name.startswith("data."):
        args.dataFileName = args.FILE
        args.dumpData = True

    elif name == "md.idx":
        args.idxFileName = args.FILE
        args.dumpIdx = True

    elif name.startswith("md."):
        args.metadataFileName = args.FILE
        args.dumpMetadata = True

    elif name.startswith("mmd."):
        args.metametadataFileName = args.FILE
        args.dumpMetaMetadata = True


def DumpIndexTableFile(args):
    global MetadataIndexTable
    global WriterMap
    global status
    indexFileList = glob.glob(args.idxFileName)
    if len(indexFileList) > 0:
        status, MetadataIndexTable, WriterMap = DumpIndexTable(
            indexFileList[0], True)
    else:
        print("There is  no BP% Index Table file as " + args.idxFileName)
        status = False
    return status


def DumpMetaMetadataFiles(args):
    global status
    mdFileList = glob.glob(args.metametadataFileName)
    if len(mdFileList) > 0:
        for fname in mdFileList:
            status = DumpMetaMetaData(fname)
    else:
        print("There are no BP% MetaMetadata files in   " +
              args.metametadataFileName)
        status = False
    return status


# xxx/md.X to xxx/md.idx
def GetIndexFileName(MDFileName):
    return MDFileName.rsplit('.', 1)[0] + ".idx"


def DumpMetadataFiles(args):
    global MetadataIndexTable
    global WriterMap
    global status
    mdFileList = glob.glob(args.metadataFileName)
    if len(mdFileList) > 0:
        if len(MetadataIndexTable) == 0:
            # need to parse index first
            IndexFileName = GetIndexFileName(mdFileList[0])
            status, MetadataIndexTable, WriterMap = DumpIndexTable(
                IndexFileName, False)

        if status:
            for fname in mdFileList:
                DumpMetaData(fname, MetadataIndexTable, WriterMap)
    else:
        print("There are no BP% Metadata files in   " + args.metadataFileName)
        status = False
    return status

# def DumpDataFiles(args):
#    dataFileList = glob.glob(args.dataFileName)
#    if len(dataFileList) > 0:
#        for fname in dataFileList:
#            DumpData(fname)
#    else:
#        print("There are no BP5 Data files in       " + args.dataFileName)


if __name__ == "__main__":

    args = SetupArgs()
    CheckFileName(args)
    # print(args)

    if args.dumpIdx:
        status = DumpIndexTableFile(args)

    if args.dumpMetaMetadata and status:
        status = DumpMetaMetadataFiles(args)

    if args.dumpMetadata and status:
        status = DumpMetadataFiles(args)

#    if args.dumpData:
#        DumpDataFiles(args)
