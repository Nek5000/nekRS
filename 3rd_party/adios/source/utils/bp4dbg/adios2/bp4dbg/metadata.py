import numpy as np
from os import fstat
from .utils import *


def ReadEncodedStringFromBuffer(buf, pos, ID, limit, lenbytes=2):
    # 'lenbytes' bytes length + string without \0
    if lenbytes == 1:
        dt = np.dtype(np.uint8)
    else:
        dt = np.dtype(np.uint16)
    namelen = np.frombuffer(buf, dtype=dt, count=1, offset=pos)[0]
    pos = pos + lenbytes
    if namelen > limit - lenbytes:
        print("ERROR: " + ID + " string length ({0}) is longer than the "
              "limit to stay inside the block ({1})".format(
                  namelen, limit - lenbytes))
        return False, "", namelen, pos
    name = buf[pos:pos + namelen].decode('ascii')
    pos = pos + namelen
    return True, name, namelen, pos


def ReadEncodedStringArrayFromBuffer(buf, pos, ID, limit, nStrings):
    s = []
    for i in range(nStrings):
        # 2 bytes length + string without \0
        namelen = np.frombuffer(buf, dtype=np.uint16, count=1, offset=pos)[0]
        pos = pos + 2
        if namelen > limit - 2:
            print("ERROR: " + ID + " string length ({0}) is longer than the "
                  "limit to stay inside the block ({1})".format(
                      namelen, limit - 2))
            return False, s, pos
        name = buf[pos:pos + namelen].decode('ascii')
        pos = pos + namelen
        limit = limit - namelen - 2
        s.append(name)
    return True, s, pos


def ReadDimensionCharacteristics(buf, pos):
    ndim = np.frombuffer(buf, dtype=np.uint8, count=1, offset=pos)[0]
    pos = pos + 1
    lgo = np.zeros(ndim, dtype=np.uint64)
    dimLen = np.frombuffer(buf, dtype=np.uint16, count=1, offset=pos)[0]
    pos = pos + 2
    if dimLen != 24 * ndim:
        print("ERROR: Encoded dimension length expected size = {0} bytes, "
              "but found {1} bytes".format(24 * ndim, dimLen))
        return False, pos, ndim, lgo

    lgo = np.frombuffer(buf, dtype=np.uint64, count=3 * ndim, offset=pos)
    pos = pos + 24 * ndim
    return True, pos, ndim, lgo


def bDataToNumpyArray(cData, typeName, nElements, startPos=0):
    if typeName == 'byte':
        return np.frombuffer(cData, dtype=np.int8, count=nElements,
                             offset=startPos)
    elif typeName == 'char':
        return np.frombuffer(cData, dtype=np.uint8, count=nElements,
                             offset=startPos)
    elif typeName == 'short':
        return np.frombuffer(cData, dtype=np.int16, count=nElements,
                             offset=startPos)
    elif typeName == 'integer':
        return np.frombuffer(cData, dtype=np.int32, count=nElements,
                             offset=startPos)
    elif typeName == 'long':
        return np.frombuffer(cData, dtype=np.int64, count=nElements,
                             offset=startPos)

    elif typeName == 'unsigned_byte':
        return np.frombuffer(cData, dtype=np.uint8, count=nElements,
                             offset=startPos)
    elif typeName == 'unsigned_short':
        return np.frombuffer(cData, dtype=np.uint16, count=nElements,
                             offset=startPos)
    elif typeName == 'unsigned_integer':
        return np.frombuffer(cData, dtype=np.uint32, count=nElements,
                             offset=startPos)
    elif typeName == 'unsigned_long':
        return np.frombuffer(cData, dtype=np.uint64, count=nElements,
                             offset=startPos)

    elif typeName == 'real':
        return np.frombuffer(cData, dtype=np.float32, count=nElements,
                             offset=startPos)
    elif typeName == 'double':
        return np.frombuffer(cData, dtype=np.float64, count=nElements,
                             offset=startPos)
    elif typeName == 'long_double':
        return np.frombuffer(cData, dtype=np.float128, count=nElements,
                             offset=startPos)

    elif typeName == 'complex':
        return np.frombuffer(cData, dtype=np.complex64, count=nElements,
                             offset=startPos)
    elif typeName == 'double_complex':
        return np.frombuffer(cData, dtype=np.complex128, count=nElements,
                             offset=startPos)

    else:
        return np.zeros(1, dtype=np.uint32)


def ReadCharacteristicsFromMetaData(buf, idx, pos, limit, typeID,
                                    fileOffset, isVarCharacteristics):
    cStartPosition = pos
    dataTypeName = GetTypeName(typeID)
    print("        Block {0}: ".format(idx))
    print("            Starting offset : {0}".format(fileOffset))
    # 1 byte NCharacteristics
    nChars = np.frombuffer(buf, dtype=np.uint8, count=1, offset=pos)[0]
    pos = pos + 1
    print("            # of Characteristics    : {0}".format(nChars))
    # 4 bytes length
    charLen = np.frombuffer(buf, dtype=np.uint8, count=32, offset=pos)[0]
    pos = pos + 4
    print("            Characteristics Length  : {0}".format(charLen))

    # For attributes, we need to remember the dimensions and size
    # when reading the value
    ndim = 0
    nElems = 1

    for i in range(nChars):
        print("            Characteristics[{0}]".format(i))
        # 1 byte TYPE
        cID = np.frombuffer(buf, dtype=np.uint8, count=1, offset=pos)[0]
        pos = pos + 1
        cName = GetCharacteristicName(cID)
        print("                Type           : {0} ({1}) ".format(
            cName, cID))
        cLen = GetCharacteristicDataLength(cID, typeID)

        if cName == 'dimensions':
            status, pos, ndim, lgo = ReadDimensionCharacteristics(buf, pos)
            if not status:
                return status, pos
            print("                # of Dims      : {0}".format(ndim))
            if ndim > 0:
                print("                Dims (lgo)     : (", end="")
                for d in range(ndim):
                    p = 3 * d
                    nElems = int(nElems * lgo[p])  # need for value later
                    print("{0}:{1}:{2}".format(lgo[p], lgo[p + 1], lgo[p + 2]),
                          end="")
                    if d < ndim - 1:
                        print(", ", end="")
                    else:
                        print(")")

        elif cName == 'value' or cName == 'min' or cName == 'max':
            if dataTypeName == 'string':
                namelimit = limit - (pos - cStartPosition)
                status, s, sLen, pos = ReadEncodedStringFromBuffer(
                    buf, pos, "String Value", namelimit)
                if not status:
                    return False, pos
                print("                Value          : '" + s +
                      "' ({0} bytes)".format(sLen))
            elif dataTypeName == 'string_array':
                namelimit = limit - (pos - cStartPosition)
                status, strList, pos = ReadEncodedStringArrayFromBuffer(
                    buf, pos, "String Array", namelimit, lgo[0])
                if not status:
                    return False, pos
                print("                Value          : [", end="")
                for j in range(len(strList)):
                    print("'" + strList[j] + "'", end="")
                    if j < len(strList) - 1:
                        print(", ", end="")
                print("]")

            else:
                if isVarCharacteristics:
                    cData = buf[pos:pos + cLen]
                    pos = pos + cLen
                    data = bDataToNumpyArray(cData, dataTypeName, 1)
                    print("                Value          : {0}"
                          "  ({1} bytes)".format(data[0], cLen))
                else:  # attribute value characteristics are different
                    dataTypeSize = GetTypeSize(typeID)
                    nBytes = int(nElems * dataTypeSize)
                    cData = buf[pos:pos + nBytes]
                    pos = pos + nBytes
                    data = bDataToNumpyArray(cData, dataTypeName, nElems)
                    print("                Value          : [", end="")
                    for j in range(nElems):
                        print("{0}".format(data[j]), end="")
                        if j < nElems - 1:
                            print(", ", end="")
                    print("]")

        elif cName == 'offset' or cName == 'payload_offset':
            cData = buf[pos:pos + cLen]
            pos = pos + cLen
            data = bDataToNumpyArray(cData, 'unsigned_long', 1)
            print("                Value          : {0}  ({1} bytes)".format(
                  data[0], cLen))
        elif cName == 'time_index' or cName == 'file_index':
            cData = buf[pos:pos + cLen]
            pos = pos + cLen
            data = bDataToNumpyArray(cData, 'unsigned_integer', 1)
            print("                Value          : {0}  ({1} bytes)".format(
                data[0], cLen))
        elif cName == 'minmax':
            nBlocks = np.frombuffer(
                buf, dtype=np.uint16, count=1, offset=pos)[0]
            print("                nBlocks        : {0}".format(nBlocks))
            pos = pos + 2
            bminmax = bDataToNumpyArray(buf, dataTypeName, 2, pos)
            pos = pos + 2 * cLen
            print("                Min/max        : {0} / {1}".format(
                bminmax[0], bminmax[1]))
            if nBlocks > 1:
                method = np.frombuffer(buf, dtype=np.uint8,
                                       count=1, offset=pos)[0]
                pos = pos + 1
                print("                Division method: {0}".format(method))
                blockSize = np.frombuffer(buf, dtype=np.uint64,
                                          count=1, offset=pos)[0]
                pos = pos + 8
                print("                Block size     : {0}".format(blockSize))
                div = np.frombuffer(buf, dtype=np.uint16,
                                    count=ndim, offset=pos)
                pos = pos + 2 * ndim
                print("                Division vector: (", end="")
                for d in range(ndim):
                    print("{0}".format(div[d]), end="")
                    if d < ndim - 1:
                        print(", ", end="")
                    else:
                        print(")")
                minmax = bDataToNumpyArray(buf, dataTypeName, 2 * nBlocks, pos)
                pos = pos + 2 * nBlocks * cLen
                for i in range(nBlocks):
                    print("                Min/max        : {0} / {1}".format(
                        minmax[2 * i], minmax[2 * i + 1]))
        elif cName == "transform_type":
            # Operator name (8 bit length)
            namelimit = limit - (pos - cStartPosition)
            status, s, sLen, pos = ReadEncodedStringFromBuffer(
                buf, pos, "Operator Name", namelimit, lenbytes=1)
            if not status:
                return False, pos
            print("                Operator       : '" + s +
                  "' ({0} bytes)".format(sLen))

            # 1 byte TYPE
            typeID = buf[pos]
            pos = pos + 1
            print("                Pre-type       : {0} ({1}) ".format(
                GetTypeName(typeID), typeID))

            # Pre-transform dimenstions
            status, pos, ndim, lgo = ReadDimensionCharacteristics(buf, pos)
            if not status:
                return status, pos
            print("                Pre-# of dims  : {0}".format(ndim))
            if ndim > 0:
                print("                Pre-Dims (lgo) : (", end="")
                for d in range(ndim):
                    p = 3 * d
                    nElems = int(nElems * lgo[p])  # need for value later
                    print("{0}:{1}:{2}".format(lgo[p], lgo[p + 1], lgo[p + 2]),
                          end="")
                    if d < ndim - 1:
                        print(", ", end="")
                    else:
                        print(")")

            # Operator specific metadata
            omdlen = np.frombuffer(buf, dtype=np.uint16,
                                   count=1, offset=pos)[0]
            pos = pos + 2
            print("                Op. data length: {0}".format(omdlen))
            pos = pos + omdlen
        else:
            print("                ERROR: could not understand this "
                  "characteristics type '{0}' id {1}".format(cName, cID))
    return True, pos


def ReadPGMD(buf, idx, pos, limit, pgStartOffset):
    # Read one PG index group
    pgStartPosition = pos
    print("    PG {0}: ".format(idx))
    print("        Starting offset : {0}".format(pgStartOffset))

    # 2 bytes PG Length + Name length
    pgLength = np.frombuffer(buf, dtype=np.uint16, count=1, offset=pos)[0]
    pos = pos + 2
    print(
        "        PG length       : {0} bytes (+2 for length)".format(
            pgLength))
    if pgStartPosition + pgLength + 2 > limit:
        print("ERROR: There is not enough bytes {0} left in PG index block "
              "to read this single PG index ({1} bytes)").format(
                  limit - pgStartPosition, pgLength)
        return False, pos

    pgNameLen = np.frombuffer(buf, dtype=np.uint16, count=1, offset=pos)[0]
    pos = pos + 2
    if pgStartPosition + pgNameLen > limit:
        print("ERROR: There is not enough bytes {0} left in PG index block "
              "to read the name of this single PG index ({1} bytes)").format(
                  limit - pos + 2, pgNameLen)
        return False, pos

    pgName = buf[pos:pos + pgNameLen].decode('ascii')
    pos = pos + pgNameLen
    print("        PG Name         : '" + pgName +
          "' ({0} bytes)".format(pgNameLen))

    # ColumnMajor (host language Fortran) 1 byte, 'y' or 'n'
    isColumnMajor = buf[pos]  # this is an integer value
    pos = pos + 1
    if isColumnMajor != ord('y') and isColumnMajor != ord('n'):
        print(
            "ERROR: Next byte for isColumnMajor must be 'y' or 'n' "
            "but it isn't = {0} (={1})".format(
                chr(isColumnMajor), isColumnMajor))
        return False
    print("        isColumnMajor   : " + chr(isColumnMajor))

    processID = np.frombuffer(buf, dtype=np.uint32, count=1, offset=pos)[0]
    pos = pos + 4
    print("        process ID      : {0}".format(processID))

    pgTimeNameLen = np.frombuffer(buf, dtype=np.uint16, count=1, offset=pos)[0]
    pos = pos + 2
    if pgStartPosition + pgTimeNameLen > limit:
        print("ERROR: There is not enough bytes {0} left in PG index block "
              "to read the name of this single PG index ({1} bytes)").format(
                  limit - pos + 2, pgTimeNameLen)
        return False, pos

    pgTimeName = buf[pos:pos + pgTimeNameLen].decode('ascii')
    pos = pos + pgTimeNameLen
    print("        Timestep Name   : '" + pgTimeName +
          "' ({0} bytes)".format(pgTimeNameLen))

    step = np.frombuffer(buf, dtype=np.uint32, count=1, offset=pos)[0]
    pos = pos + 4
    print("        Step            : {0}".format(step))

    ptr = np.frombuffer(buf, dtype=np.uint64, count=1, offset=pos)[0]
    pos = pos + 8
    print("        Offset in data  : {0}".format(ptr))

    return True, pos


def ReadVarMD(buf, idx, pos, limit, varStartOffset):
    # Read one VAR index group
    varStartPosition = pos
    print("    Var {0}: ".format(idx))
    print("        Starting offset : {0}".format(varStartOffset))

    # 4 bytes VAR Index Length
    varLength = np.frombuffer(buf, dtype=np.uint32, count=1, offset=pos)[0]
    pos = pos + 4
    print("        Var idx length  : {0} bytes (+4 for idx length)".format(
        varLength))
    if varStartPosition + varLength + 4 > limit:
        print("ERROR: There is not enough bytes in Var index block "
              "to read this single Var index")
        return False, pos

    memberID = np.frombuffer(buf, dtype=np.uint32, count=1, offset=pos)[0]
    pos = pos + 4
    print("        MemberID        : {0}".format(memberID))

    namelimit = limit - (pos - varStartPosition)
    status, grpName, grpNameLen, pos = ReadEncodedStringFromBuffer(
        buf, pos, "Group Name", namelimit)
    if not status:
        return False, pos
    print("        Group Name      : '" + grpName +
          "' ({0} bytes)".format(grpNameLen))

    namelimit = limit - (pos - varStartPosition)
    status, varName, varNameLen, pos = ReadEncodedStringFromBuffer(
        buf, pos, "Var Name", namelimit)
    if not status:
        return False, pos
    print("        Var Name        : '" + varName +
          "' ({0} bytes)".format(varNameLen))

    # print("        Current offset : {0}".format(varStartOffset + pos))
    # namelimit = limit - (pos - varStartPosition)
    # status, varPath, varPathLen, pos = ReadEncodedStringFromBuffer(
    #     buf, pos, "Var Path", namelimit)
    # if not status:
    #     return False, pos
    # print("        Var Path        : '" + varPath +
    #       "' ({0} bytes)".format(varPathLen))

    # 1 byte ORDER (K, C, F)
    order = buf[pos]  # this is an integer value
    pos = pos + 1
    if order != ord('K') and order != ord('C') and order != ord('F'):
        print(
            "ERROR: Next byte for Order must be 'K', 'C', or 'F' "
            "but it isn't = {0} (={1})".format(
                chr(order), order))
        return False
    print("        Order           : " + chr(order))

    # 1 byte UNUSED
    unused = buf[pos]  # this is an integer value
    pos = pos + 1
    print("        Unused byte     : {0}".format(unused))

    # 1 byte TYPE
    typeID = buf[pos]
    pos = pos + 1
    print("        Type            : {0} ({1}) ".format(
        GetTypeName(typeID), typeID))

    # 8 byte Number of Characteristics Sets
    cSets = np.frombuffer(buf, dtype=np.uint64, count=1, offset=pos)[0]
    pos = pos + 8
    print("        # of blocks     : {0}".format(cSets))

#   This loop only reads the number of reported blocks
#     for i in range(cSets):
#         # one characteristics block
#         newlimit = limit - (pos - varStartPosition)
#         fileOffset = varStartOffset + (pos - varStartPosition)
#         status, pos = ReadCharacteristicsFromMetaData(
#             buf, i, pos, newlimit, typeID, fileOffset, True)
#         if not status:
#             return False

#   This loop reads blocks until the reported length of variable index length
    i = 0
    while pos < varStartPosition + varLength:
        # one characteristics block
        newlimit = limit - (pos - varStartPosition)
        fileOffset = varStartOffset + (pos - varStartPosition)
        status, pos = ReadCharacteristicsFromMetaData(
            buf, i, pos, newlimit, typeID, fileOffset, True)
        if not status:
            return False
        i = i + 1

    if (i != cSets):
        print(
            "ERROR: reported # of blocks (={0}) != # of encoded blocks "
            " (={1})".format(
                cSets, i))

    return True, pos


def ReadAttrMD(buf, idx, pos, limit, attrStartOffset):
    # Read one ATTR index group
    attrStartPosition = pos
    print("    Attr {0}: ".format(idx))
    print("        Starting offset : {0}".format(attrStartOffset))

    # 4 bytes ATTR Index Length
    attrLength = np.frombuffer(buf, dtype=np.uint32, count=1, offset=pos)[0]
    pos = pos + 4
    print("        Attr idx length : {0} bytes (+4 for idx length)".format(
        attrLength))
    if attrStartPosition + attrLength + 4 > limit:
        print("ERROR: There is not enough bytes in Attr index block "
              "to read this single Attr index")
        return False, pos

    memberID = np.frombuffer(buf, dtype=np.uint32, count=1, offset=pos)[0]
    pos = pos + 4
    print("        MemberID        : {0}".format(memberID))

    namelimit = limit - (pos - attrStartPosition)
    status, grpName, grpNameLen, pos = ReadEncodedStringFromBuffer(
        buf, pos, "Group Name", namelimit)
    if not status:
        return False, pos
    print("        Group Name      : '" + grpName +
          "' ({0} bytes)".format(grpNameLen))

    namelimit = limit - (pos - attrStartPosition)
    status, attrName, attrNameLen, pos = ReadEncodedStringFromBuffer(
        buf, pos, "Attr Name", namelimit)
    if not status:
        return False, pos
    print("        Attr Name       : '" + attrName +
          "' ({0} bytes)".format(attrNameLen))

    # print("        Current offset : {0}".format(attrStartOffset + pos))
    namelimit = limit - (pos - attrStartPosition)
    status, attrPath, attrPathLen, pos = ReadEncodedStringFromBuffer(
        buf, pos, "Attr Path", namelimit)
    if not status:
        return False, pos
    print("        Attr Path       : '" + attrPath +
          "' ({0} bytes)".format(attrPathLen))

    # 1 byte TYPE
    typeID = buf[pos]
    pos = pos + 1
    print("        Type            : {0} ({1}) ".format(
        GetTypeName(typeID), typeID))

    # 8 byte Number of Characteristics Sets
    cSets = np.frombuffer(buf, dtype=np.uint64, count=1, offset=pos)[0]
    pos = pos + 8
    print("        # of blocks     : {0}".format(cSets))

    for i in range(cSets):
        # one characteristics block
        newlimit = limit - (pos - attrStartPosition)
        fileOffset = attrStartOffset + (pos - attrStartPosition)
        status, pos = ReadCharacteristicsFromMetaData(
            buf, i, pos, newlimit, typeID, fileOffset, False)
        if not status:
            return False

    return True, pos


def ReadMetadataStep(f, fileSize, step):
    # Read metadata of one step
    mdStartPosition = f.tell()
    if step > 0:
        print("========================================================")
    print("Step {0}: ".format(step))
    print("    PG Index offset   : {0}".format(mdStartPosition))

    # Read the PG Index

    # 8 bytes PG Count + Index Length
    pgInfo = np.fromfile(f, dtype=np.uint64, count=2)
    pgCount = pgInfo[0]
    pgLength = pgInfo[1]
    print("    # of PGs          : {0}".format(pgCount))
    print("    PG Index length   : {0} bytes (+16 for count+len)".format(
        pgLength))

    pgStartPosition = f.tell()
    if pgStartPosition + pgLength > fileSize:
        print("ERROR: There is not enough bytes in file to read the PG index")
        return False

    pgmd = f.read(pgLength)
    pgmdPos = 0
    for i in range(pgCount):
        # VMD block
        status, pgmdPos = ReadPGMD(
            pgmd, i, pgmdPos, pgLength, pgStartPosition + pgmdPos)
        if not status:
            return False

    # Read the VAR Index

    print("    ----------------------------------")
    print("    Var Index offset  : {0}".format(f.tell()))

    # 4 bytes VAR Count + 8 bytes Var Index Length
    varCount = np.fromfile(f, dtype=np.uint32, count=1)[0]
    varLength = np.fromfile(f, dtype=np.uint64, count=1)[0]
    print("    # of Variables    : {0}".format(varCount))
    print("    Var Index length  : {0} bytes (+12 for count+len)".format(
        varLength))

    varsStartPosition = f.tell()
    if varsStartPosition + varLength > fileSize:
        print("ERROR: There is not enough bytes in file to read the VAR index")
        return False

    varmd = f.read(varLength)
    varmdPos = 0
    for i in range(varCount):
        # VMD block
        status, varmdPos = ReadVarMD(
            varmd, i, varmdPos, varLength, varsStartPosition + varmdPos)
        if not status:
            return False

    # Read the ATTR Index

    print("    ----------------------------------")
    print("    Attr Index offset : {0}".format(f.tell()))

    # 4 bytes ATTR Count + 8 bytes Attr Index Length
    attrCount = np.fromfile(f, dtype=np.uint32, count=1)[0]
    attrLength = np.fromfile(f, dtype=np.uint64, count=1)[0]
    print("    # of attriables    : {0}".format(attrCount))
    print("    Attr Index length  : {0} bytes (+12 for count+len)".format(
        attrLength))

    attrsStartPosition = f.tell()
    if attrsStartPosition + attrLength > fileSize:
        print("ERROR: There is not enough bytes in file "
              "to read the Attribute index")
        return False

    attrmd = f.read(attrLength)
    attrmdPos = 0
    for i in range(attrCount):
        # VMD block
        status, attrmdPos = ReadAttrMD(
            attrmd, i, attrmdPos, attrLength, attrsStartPosition + attrmdPos)
        if not status:
            return False

    return True


def DumpMetaData(fileName):
    print("========================================================")
    print("    Metadata File: " + fileName)
    print("========================================================")
    with open(fileName, "rb") as f:
        fileSize = fstat(f.fileno()).st_size
        status = ReadHeader(f, fileSize, "Metadata")
        step = 0
        while (f.tell() < fileSize - 12 and status):
            status = ReadMetadataStep(f, fileSize, step)
            step = step + 1
    return status


if __name__ == "__main__":
    print("ERROR: Utility main program is bp4dbg.py")
