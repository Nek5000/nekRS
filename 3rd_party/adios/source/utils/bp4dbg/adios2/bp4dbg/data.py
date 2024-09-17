import numpy as np
from os import fstat
from .utils import *


def ReadEncodedString(f, ID, limit, lensize=2):
    if lensize == 2:
        # 2 bytes length + string without \0
        namelen = np.fromfile(f, dtype=np.uint16, count=1)[0]
    elif lensize == 4:
        # 2 bytes length + string without \0
        namelen = np.fromfile(f, dtype=np.uint32, count=1)[0]
    else:
        print("CODING ERROR: bp4dbp_data.ReadEncodedString: "
              "lensize must be 2 or 4")
        return False, ""
    if namelen > limit:
        print("ERROR: " + ID + " string length ({0}) is longer than the "
              "limit to stay inside the block ({1})".format(
                  namelen, limit))
        return False, ""
    name = f.read(namelen).decode('ascii')
    return True, name


def ReadEncodedStringArray(f, ID, limit, nStrings):
    s = []
    for i in range(nStrings):
        # 2 bytes length + string
        # !!! String here INCLUDES Terminating \0 !!!
        namelen = np.fromfile(f, dtype=np.uint32, count=1)[0]
        if namelen > limit - 4:
            print("ERROR: " + ID + " string length ({0}) is longer than the "
                  "limit to stay inside the block ({1})".format(
                      namelen, limit - 4))
            return False, s
        name = f.read(namelen).decode('ascii')
        limit = limit - namelen - 4
        s.append(name[0:-1])  # omit the terminating \0
    return True, s


def readDataToNumpyArray(f, typeName, nElements):
    if typeName == 'byte':
        return np.fromfile(f, dtype=np.int8, count=nElements)
    elif typeName == 'char':
        return np.fromfile(f, dtype=np.uint8, count=nElements)
    elif typeName == 'short':
        return np.fromfile(f, dtype=np.int16, count=nElements)
    elif typeName == 'integer':
        return np.fromfile(f, dtype=np.int32, count=nElements)
    elif typeName == 'long':
        return np.fromfile(f, dtype=np.int64, count=nElements)

    elif typeName == 'unsigned_byte':
        return np.fromfile(f, dtype=np.uint8, count=nElements)
    elif typeName == 'unsigned_short':
        return np.fromfile(f, dtype=np.uint16, count=nElements)
    elif typeName == 'unsigned_integer':
        return np.fromfile(f, dtype=np.uint32, count=nElements)
    elif typeName == 'unsigned_long':
        return np.fromfile(f, dtype=np.uint64, count=nElements)

    elif typeName == 'real':
        return np.fromfile(f, dtype=np.float32, count=nElements)
    elif typeName == 'double':
        return np.fromfile(f, dtype=np.float64, count=nElements)
    elif typeName == 'long_double':
        return np.fromfile(f, dtype=np.float128, count=nElements)

    elif typeName == 'complex':
        return np.fromfile(f, dtype=np.complex64, count=nElements)
    elif typeName == 'double_complex':
        return np.fromfile(f, dtype=np.complex128, count=nElements)

    else:
        return np.zeros(1, dtype=np.uint32)


def ReadCharacteristicsFromData(f, limit, typeID, ndim):
    cStartPosition = f.tell()
    dataTypeName = GetTypeName(typeID)
    # 1 byte NCharacteristics
    nCharacteristics = np.fromfile(f, dtype=np.uint8, count=1)[0]
    print("      # of Characteristics    : {0}".format(nCharacteristics))
    # 4 bytes length
    charLen = np.fromfile(f, dtype=np.uint32, count=1)[0]
    print("      Characteristics Length  : {0}".format(charLen))

    for i in range(nCharacteristics):
        print("      Characteristics[{0}]".format(i))
        # 1 byte TYPE
        cID = np.fromfile(f, dtype=np.uint8, count=1)[0]
        cName = GetCharacteristicName(cID)
        print("          Type        : {0} ({1}) ".format(cName, cID))
        if cName == 'value' or cName == 'min' or cName == 'max':
            if dataTypeName == 'string':
                namelimit = limit - (f.tell() - cStartPosition)
                status, s = ReadEncodedString(f, "String Value", namelimit)
                if not status:
                    return False
                print("          Value       : '" + s + "'")
            else:
                data = readDataToNumpyArray(f, dataTypeName, 1)
                print("          Value       : {0}".format(data[0]))
        elif cName == 'offset' or cName == 'payload_offset':
            data = readDataToNumpyArray(f, 'unsigned_long', 1)
            print("          Value       : {0}".format(data[0]))
        elif cName == 'time_index' or cName == 'file_index':
            data = readDataToNumpyArray(f, 'unsigned_integer', 1)
            print("          Value       : {0}".format(data[0]))
        elif cName == 'minmax':
            nBlocks = np.fromfile(f,
                                  dtype=np.uint16, count=1)[0]
            print("          nBlocks     : {0}".format(nBlocks))
            bminmax = readDataToNumpyArray(f, dataTypeName, 2)
            print("          Min/max     : {0} / {1}".format(
                bminmax[0], bminmax[1]))
            if nBlocks > 1:
                method = np.fromfile(f, dtype=np.uint8,
                                     count=1)[0]
                print("          Division method: {0}".format(method))
                blockSize = np.fromfile(f, dtype=np.uint64,
                                        count=1)[0]
                print("          Block size     : {0}".format(blockSize))
                div = np.fromfile(f, dtype=np.uint16,
                                  count=ndim)
                print("          Division vector: (", end="")
                for d in range(ndim):
                    print("{0}".format(div[d]), end="")
                    if d < ndim - 1:
                        print(", ", end="")
                    else:
                        print(")")
                minmax = readDataToNumpyArray(
                    f, dataTypeName, 2 * nBlocks)
                for i in range(nBlocks):
                    print("          Min/max        : {0} / {1}".format(
                        minmax[2 * i], minmax[2 * i + 1]))
        else:
            print("                ERROR: could not understand this "
                  "characteristics type '{0}' id {1}".format(cName, cID))
    return True

# Read String Variable data


def ReadStringVarData(f, expectedSize,
                      varsStartPosition):
    # 2 bytes String Length
    len = np.fromfile(f, dtype=np.uint16, count=1)[0]
    if len != expectedSize - 2:
        print("ERROR: Variable data block size does not equal the size "
              "calculated from var block length")
        print("Expected size = {0}  calculated size "
              "from encoded length info {1}".
              format(expectedSize, len + 2))
        return False

    str = f.read(len).decode('ascii')
    print("      Variable Data   : '" + str + "'")
    return True

# Read Variable data


def ReadVarData(f, nElements, typeID, ldims, varLen,
                varsStartPosition, varsTotalLength):
    if typeID == 9:  # string type
        return ReadStringVarData(f, varLen, varsStartPosition)
    typeSize = GetTypeSize(typeID)
    if typeSize == 0:
        print("ERROR: Cannot process variable data block with "
              "unknown type size")
        return False

    currentPosition = f.tell()
    print("      Payload offset  : {0}".format(currentPosition))

    if currentPosition + varLen > varsStartPosition + varsTotalLength:
        print("ERROR: Variable data block of size would reach beyond all "
              "variable blocks")
        print("VarsStartPosition = {0} varsTotalLength = {1}".format(
            varsStartPosition, varsTotalLength))
        print("current Position = {0} var block length = {1}".format(
            currentPosition, varLen))
        return False

    nBytes = int(varLen.item())

    if nElements == 1:
        # single value. read and print
        value = readDataToNumpyArray(f, GetTypeName(typeID),
                                     nElements)
        print("      Payload (value) : {0} ({1} bytes)".format(
            value[0], nBytes))
    else:
        # seek instead of reading for now
        # f.read(nBytes)
        f.seek(nBytes, 1)
        # data = readDataToNumpyArray(f, GetTypeName(typeID),
        #                            nElements)
        print("      Payload (array) : {0} bytes".format(nBytes))

    return True

# Read a variable's metadata


def ReadVMD(f, varidx, varsStartPosition, varsTotalLength):
    startPosition = f.tell()
    print("  Var {0:5d}".format(varidx))
    print("      Starting offset : {0}".format(startPosition))
    # 4 bytes TAG
    tag = f.read(4)
    if tag != b"[VMD":
        print("  Tag: " + str(tag))
        print("ERROR: VAR group does not start with [VMD")
        return False
    print("      Tag             : " + tag.decode('ascii'))

    # 8 bytes VMD Length
    vmdlen = np.fromfile(f, dtype=np.uint64, count=1)[0]
    print("      Var block size  : {0} bytes (+4 for Tag)".format(vmdlen))
    expectedVarBlockLength = vmdlen + 4  # [VMD is not included in vmdlen

    if startPosition + expectedVarBlockLength > \
            varsStartPosition + varsTotalLength:
        print("ERROR: There is not enough bytes inside this PG to read "
              "this Var block")
        print("VarsStartPosition = {0} varsTotalLength = {1}".format(
            varsStartPosition, varsTotalLength))
        print("current var's start position = {0} var block length = {1}".
              format(startPosition, expectedVarBlockLength))
        return False

    # 4 bytes VAR MEMBER ID
    memberID = np.fromfile(f, dtype=np.uint32, count=1)[0]
    print("      Member ID       : {0}".format(memberID))

    # VAR NAME, 2 bytes length + string without \0
    sizeLimit = expectedVarBlockLength - (f.tell() - startPosition)
    status, varname = ReadEncodedString(f, "Var Name", sizeLimit)
    if not status:
        return False
    print("      Var Name        : " + varname)

    # VAR PATH, 2 bytes length + string without \0
    # sizeLimit = expectedVarBlockLength - (f.tell() - startPosition)
    # status, varpath = ReadEncodedString(f, "Var Path", sizeLimit)
    # if not status:
    #     return False
    # print("      Var Path        : " + varpath)

    # 1 byte ORDER (K, C, F)
    order = f.read(1)
    if (order != b'K' and order != b'C' and order != b'F' and order != b'\x00'):
        print(
            "ERROR: Next byte for Order must be 'K', 'C', or 'F' "
            "but it isn't = {0}".format(order))
        return False
    if (order == b'\x00'):
        order = b'0'
    print("        Order           : " + order.decode('ascii'))

    # 1 byte UNUSED
    unused = f.read(1)
    print("        Unused byte     : {0}".format(ord(unused)))

    # 1 byte TYPE
    typeID = np.fromfile(f, dtype=np.uint8, count=1)[0]
    print("      Type            : {0} ({1}) ".format(
        GetTypeName(typeID), typeID))

    # ISDIMENSIONS 1 byte, 'y' or 'n'
    isDimensionVar = f.read(1)
    if (isDimensionVar != b'y' and isDimensionVar != b'n'):
        print(
            "ERROR: Next byte for isDimensionVar must be 'y' or 'n' "
            "but it isn't = {0}".format(isDimensionVar))
        return False
    print("      isDimensionVar  : " + isDimensionVar.decode('ascii'))

    # 1 byte NDIMENSIONS
    ndims = np.fromfile(f, dtype=np.uint8, count=1)[0]
    print("      # of Dimensions : {0}".format(
        ndims))

    # DIMLENGTH
    dimsLen = np.fromfile(f, dtype=np.uint16, count=1)[0]
    print("      Dims Length     : {0}".format(
        dimsLen))

    nElements = np.uint64(1)
    ldims = np.zeros(ndims, dtype=np.uint64)
    isLocalValueArray = False
    for i in range(ndims):
        print("      Dim[{0}]".format(i))
        # Read Local Dimensions (1 byte flag + 8 byte value)
        # Is Dimension a variable ID 1 byte, 'y' or 'n' or '\0'
        isDimensionVarID = f.read(1)
        if isDimensionVarID != b'y' and isDimensionVarID != b'n' and \
                isDimensionVarID != b'\0':
            print(
                "ERROR: Next byte for isDimensionVarID must be 'y' or 'n' "
                "but it isn't = {0}".format(isDimensionVarID))
            return False
        if isDimensionVarID == b'\0':
            isDimensionVarID = b'n'
        ldims[i] = np.fromfile(f, dtype=np.uint64, count=1)[0]
        print("           local  dim : {0}".format(ldims[i]))
        nElements = nElements * ldims[i]
        # Read Global Dimensions (1 byte flag + 8 byte value)
        # Is Dimension a variable ID 1 byte, 'y' or 'n' or '\0'
        isDimensionVarID = f.read(1)
        if isDimensionVarID != b'y' and isDimensionVarID != b'n' \
                and isDimensionVarID != b'\0':
            print(
                "ERROR: Next byte for isDimensionVarID must be 'y' or 'n' "
                "but it isn't = {0}".format(isDimensionVarID))
            return False
        if isDimensionVarID == b'\0':
            isDimensionVarID = b'n'
        gdim = np.fromfile(f, dtype=np.uint64, count=1)[0]
        if i == 0 and ldims[i] == 0 and gdim == LocalValueDim:
            print("           global dim : LocalValueDim ({0})".format(gdim))
            isLocalValueArray = True
        else:
            print("           global dim : {0}".format(gdim))

        # Read Offset Dimensions (1 byte flag + 8 byte value)
        # Is Dimension a variable ID 1 byte, 'y' or 'n' or '\0'
        isDimensionVarID = f.read(1)
        if isDimensionVarID != b'y' and isDimensionVarID != b'n' and \
                isDimensionVarID != b'\0':
            print(
                "ERROR: Next byte for isDimensionVarID must be 'y' or 'n' "
                "but it isn't = {0}".format(isDimensionVarID))
            return False
        if isDimensionVarID == b'\0':
            isDimensionVarID = b'n'
        offset = np.fromfile(f, dtype=np.uint64, count=1)[0]
        print("           offset dim : {0}".format(offset))

    sizeLimit = expectedVarBlockLength - (f.tell() - startPosition)
    status = ReadCharacteristicsFromData(f, sizeLimit, typeID, ndims)
    if not status:
        return False

    # Padded end TAG
    # 1 byte length of tag
    endTagLen = np.fromfile(f, dtype=np.uint8, count=1)[0]
    tag = f.read(endTagLen)
    if not tag.endswith(b"VMD]"):
        print("  Tag: " + str(tag))
        print("ERROR: VAR group metadata does not end with VMD]")
        return False
    print("      Tag (pad {0:2d})    : {1}".format(
        endTagLen - 4, tag.decode('ascii')))

    # special case: LocalValueDim: local values turned into 1D global array
    # but it seems there is no data block at all for these variables
    if isLocalValueArray:
        ldims[0] = 1
        nElements = np.uint64(1)
    else:
        expectedVarDataSize = expectedVarBlockLength - \
            (f.tell() - startPosition)
        status = ReadVarData(f, nElements, typeID, ldims, expectedVarDataSize,
                             varsStartPosition, varsTotalLength)
    if not status:
        return False

    return True

# Read an attribute's metadata and value


def ReadAMD(f, attridx, attrsStartPosition, attrsTotalLength):
    startPosition = f.tell()
    print("  attr {0:5d}".format(attridx))
    print("      Starting offset : {0}".format(startPosition))
    # 4 bytes TAG
    tag = f.read(4)
    if tag != b"[AMD":
        print("  Tag: " + str(tag))
        print("ERROR: ATTR group does not start with [AMD")
        return False
    print("      Tag             : " + tag.decode('ascii'))

    # 8 bytes AMD Length
    amdlen = np.fromfile(f, dtype=np.uint32, count=1)[0]
    print("      Attr block size : {0} bytes (+4 for Tag)".format(amdlen))
    expectedAttrBlockLength = amdlen + 4  # [AMD is not included in amdlen
    if startPosition + expectedAttrBlockLength > \
            attrsStartPosition + attrsTotalLength:
        print("ERROR: There is not enough bytes inside this PG "
              "to read this Attr block")
        print("AttrsStartPosition = {0} attrsTotalLength = {1}".format(
            attrsStartPosition, attrsTotalLength))
        print("current attr's start position = {0} "
              "attr block length = {1}".format(
                  startPosition, expectedAttrBlockLength))
        return False

    # 4 bytes ATTR MEMBER ID
    memberID = np.fromfile(f, dtype=np.uint32, count=1)[0]
    print("      Member ID       : {0}".format(memberID))

    # ATTR NAME, 2 bytes length + string without \0
    sizeLimit = expectedAttrBlockLength - (f.tell() - startPosition)
    status, attrname = ReadEncodedString(f, "Attr Name", sizeLimit)
    if not status:
        return False
    print("      Attr Name       : " + attrname)

    # ATTR PATH, 2 bytes length + string without \0
    sizeLimit = expectedAttrBlockLength - (f.tell() - startPosition)
    status, attrpath = ReadEncodedString(f, "Attr Path", sizeLimit)
    if not status:
        return False
    print("      Attr Path       : " + attrpath)

    # isAttrAVar 1 byte, 'y' or 'n'
    isAttrAVar = f.read(1)
    if isAttrAVar != b'y' and isAttrAVar != b'n':
        print(
            "ERROR: Next byte for isAttrAVar must be 'y' or 'n' "
            "but it isn't = {0}".format(isAttrAVar))
        return False
    print("      Refers to Var?  : " + isAttrAVar.decode('ascii'))

    # 1 byte TYPE
    typeID = np.fromfile(f, dtype=np.uint8, count=1)[0]
    typeName = GetTypeName(typeID)
    print("      Type            : {0} ({1}) ".format(typeName, typeID))

    # Read Attribute data
    if typeName == 'string':
        sizeLimit = expectedAttrBlockLength - (f.tell() - startPosition)
        status, s = ReadEncodedString(
            f, "Attribute String Value", sizeLimit, 4)
        if not status:
            return False
        print("      Value           : '" + s + "'")

    elif typeName == 'string_array':
        nElems = np.fromfile(f, dtype=np.uint32, count=1)[0]
        sizeLimit = expectedAttrBlockLength - (f.tell() - startPosition)
        status, strList = ReadEncodedStringArray(
            f, "Attribute String Array", sizeLimit, nElems)
        if not status:
            return False
        print("      Value           : [", end="")
        for j in range(len(strList)):
            print("'" + strList[j] + "'", end="")
            if j < len(strList) - 1:
                print(", ", end="")
        print("]")
    else:
        nBytes = np.fromfile(f, dtype=np.uint32, count=1)[0]
        typeSize = GetTypeSize(typeID)
        nElems = int(nBytes / typeSize)
        data = readDataToNumpyArray(f, typeName, nElems)
        print("      Value           : [", end="")
        for j in range(nElems):
            print("{0}".format(data[j]), end="")
            if j < nElems - 1:
                print(", ", end="")
        print("]")

    # End TAG AMD]
    tag = f.read(4)
    if tag != b"AMD]":
        print("  Tag: " + str(tag))
        print("ERROR: PG group metadata does not end with AMD]")
        return False
    print("      Tag             : {0}".format(tag.decode('ascii')))

    return True

# Read one PG process group (variables and attributes from one process in
# one step)


def ReadPG(f, fileSize, pgidx):
    pgStartPosition = f.tell()
    if pgidx > 0:
        print("========================================================")
    print("Process Group {0}: ".format(pgidx))
    print("  Starting offset : {0}".format(pgStartPosition))
    tag = f.read(4)
    if tag != b"[PGI":
        print("  Tag: " + str(tag))
        print("ERROR: PG group does not start with [PGI")
        return False

    print("  Tag             : " + tag.decode('ascii'))

    # 8 bytes PG Length
    pglen = np.fromfile(f, dtype=np.uint64, count=1)[0]
    print("  PG length       : {0} bytes (+4 for Tag)".format(pglen))
    # pglen does not include the opening tag 4 bytes:
    expectedPGLength = pglen + 4
    if pgStartPosition + expectedPGLength > fileSize:
        print("ERROR: There is not enough bytes in file to read this PG")
        return False

    # ColumnMajor (host language Fortran) 1 byte, 'y' or 'n'
    isColumnMajor = f.read(1)
    if isColumnMajor != b'y' and isColumnMajor != b'n':
        print(
            "ERROR: Next byte for isColumnMajor must be 'y' or 'n' "
            "but it isn't = {0}".format(isColumnMajor))
        return False
    print("  isColumnMajor   : " + isColumnMajor.decode('ascii'))

    # PG Name, 2 bytes length + string without \0
    sizeLimit = expectedPGLength - (f.tell() - pgStartPosition)
    status, pgname = ReadEncodedString(f, "PG Name", sizeLimit)
    if not status:
        return False
    print("  PG Name         : " + pgname)

    # 4 bytes unused (for Coordination variable)
    tag = f.read(4)
    print("  Unused 4 bytes  : " + str(tag))

    # Timestep name
    sizeLimit = expectedPGLength - (f.tell() - pgStartPosition)
    status, tsname = ReadEncodedString(f, "Timestep Name", sizeLimit)
    if not status:
        return False
    print("  Step Name       : " + tsname)

    # STEP 4 bytes
    step = np.fromfile(f, dtype=np.uint32, count=1)[0]
    print("  Step Value      : {0}".format(step))

    # Method Count 1 byte1
    nMethods = np.fromfile(f, dtype=np.uint8, count=1)[0]
    print("  Methods count   : {0}".format(nMethods))

    # Method Length 2 byte1
    lenMethods = np.fromfile(f, dtype=np.uint16, count=1)[0]
    print("  Methods length  : {0}".format(lenMethods))

    print("  Methods info")
    for i in range(nMethods):
        # Method ID
        methodID = np.fromfile(f, dtype=np.uint8, count=1)[0]
        print("      Method ID   : {0}".format(methodID))
        sizeLimit = expectedPGLength - (f.tell() - pgStartPosition)
        status, methodParams = ReadEncodedString(
            f, "Method Parameters", sizeLimit)
        if not status:
            return False
        print('      M. params   : "' + methodParams + '"')

    # VARIABLES

    # VARS COUNT 4 bytes
    nVars = np.fromfile(f, dtype=np.uint32, count=1)[0]
    print("  # of Variables  : {0}".format(nVars))

    # VARS SIZE 8 bytes
    varlen = np.fromfile(f, dtype=np.uint64, count=1)[0]
    print("  Vars length     : {0} bytes".format(varlen))
    sizeLimit = expectedPGLength - (f.tell() - pgStartPosition)
    expectedVarsLength = varlen  # need to read this more
    if expectedVarsLength > sizeLimit:
        print("ERROR: There is not enough bytes in PG to read the variables")
        return False

    varsStartPosition = f.tell()
    for i in range(nVars):
        # VMD block
        status = ReadVMD(f, i, varsStartPosition, expectedVarsLength)
        if not status:
            return False

    # ATTRIBUTES

    # ATTRS COUNT 4 bytes
    nAttrs = np.fromfile(f, dtype=np.uint32, count=1)[0]
    print("  # of Attributes : {0}".format(nAttrs))

    attrsStartPosition = f.tell()
    # ATTS SIZE 8 bytes
    # attlen includes the 8 bytes of itself, so remember position before this
    attlen = np.fromfile(f, dtype=np.uint64, count=1)[0]
    print("  Attrs length    : {0} bytes".format(attlen))
    sizeLimit = expectedPGLength - (attrsStartPosition - pgStartPosition) - 4
    expectedAttrsLength = attlen  # need to read this more before reaching PGI]

    if expectedAttrsLength > sizeLimit:
        print("ERROR: There is not enough bytes in PG to read the attributes")
        return False

    attrsStartPosition = f.tell()
    for i in range(nAttrs):
        # AMD block
        status = ReadAMD(f, i, attrsStartPosition, expectedAttrsLength)
        if not status:
            return False

    # End TAG PGI]
    tag = f.read(4)
    if tag != b"PGI]":
        print("  Tag: " + str(tag))
        print("ERROR: PG group metadata does not end with PGI]")
        return False
    print("  Tag               : {0}".format(tag.decode('ascii')))

    return True


def DumpData(fileName):
    print("========================================================")
    print("    Data File: " + fileName)
    print("========================================================")
    with open(fileName, "rb") as f:
        fileSize = fstat(f.fileno()).st_size
        status = ReadHeader(f, fileSize, "Data")
        if not status:
            return status
        pgidx = 0
        while (f.tell() < fileSize - 12 and status):
            status = ReadPG(f, fileSize, pgidx)
            pgidx = pgidx + 1
    return status


if __name__ == "__main__":
    print("ERROR: Utility main program is bp4dbg.py")
