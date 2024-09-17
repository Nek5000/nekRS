/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressNull.cpp
 *
 *  Created on: Dec 1, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "CompressNull.h"
#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace core
{
namespace compress
{

CompressNull::CompressNull(const Params &parameters)
: Operator("null", COMPRESS_NULL, "compress", parameters)
{
}

size_t CompressNull::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                             const DataType varType, char *bufferOut)
{
    const uint8_t bufferVersion = 1;
    size_t bufferOutOffset = 0;
    MakeCommonHeader(bufferOut, bufferOutOffset, bufferVersion);
    size_t totalInputBytes = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(varType));
    PutParameter(bufferOut, bufferOutOffset, totalInputBytes);
    std::memcpy(bufferOut + bufferOutOffset, dataIn, totalInputBytes);
    bufferOutOffset += totalInputBytes;
    return bufferOutOffset;
}

size_t CompressNull::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    size_t bufferInOffset = 4; // skip common header

    const size_t totalBytes = GetParameter<size_t>(bufferIn, bufferInOffset);
    std::memcpy(dataOut, bufferIn + bufferInOffset, totalBytes);
    return totalBytes;
}

bool CompressNull::IsDataTypeValid(const DataType type) const { return true; }

} // end namespace compress
} // end namespace core
} // end namespace adios2
