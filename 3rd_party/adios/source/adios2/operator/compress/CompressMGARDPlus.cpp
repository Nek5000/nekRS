/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressMGARDPlus.cpp :
 *
 *  Created on: Feb 23, 2022
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "CompressMGARDPlus.h"
#include "CompressMGARD.h"
#include "adios2/core/Engine.h"
#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace core
{
namespace compress
{

CompressMGARDPlus::CompressMGARDPlus(const Params &parameters)
: Operator("mgardplus", COMPRESS_MGARDPLUS, "compress", parameters)
{
}

size_t CompressMGARDPlus::Operate(const char *dataIn, const Dims &blockStart,
                                  const Dims &blockCount, const DataType type, char *bufferOut)
{

    // Read ADIOS2 files from here
    if (!m_Parameters["MeshFile"].empty())
    {
        adios2::core::ADIOS adios("C++");
        auto &io = adios.DeclareIO("SubIO");
        auto *engine = &io.Open(m_Parameters["MeshFile"], adios2::Mode::Read);
        auto var = io.InquireVariable<float>(m_Parameters["MeshVariable"]);
        std::vector<float> data(std::accumulate(var->m_Shape.begin(), var->m_Shape.end(),
                                                sizeof(float), std::multiplies<size_t>()));
        engine->Get(*var, data);
    }
    // Read ADIOS2 files end, use data for your algorithm

    size_t bufferOutOffset = 0;
    const uint8_t bufferVersion = 1;

    MakeCommonHeader(bufferOut, bufferOutOffset, bufferVersion);

    bufferOutOffset += 32; // TODO: reserve memory space

    CompressMGARD mgard(m_Parameters);
    size_t mgardBufferSize =
        mgard.Operate(dataIn, blockStart, blockCount, type, bufferOut + bufferOutOffset);
    if (mgardBufferSize == 0)
    {
        headerSize += (bufferOutOffset + mgard.GetHeaderSize());
        return 0;
    }

    if (*reinterpret_cast<OperatorType *>(bufferOut + bufferOutOffset) == COMPRESS_MGARD)
    {
        std::vector<char> tmpDecompressBuffer(
            helper::GetTotalSize(blockCount, helper::GetDataTypeSize(type)));

        mgard.InverseOperate(bufferOut + bufferOutOffset, mgardBufferSize,
                             tmpDecompressBuffer.data());

        // TODO: now the original data is in dataIn, the compressed and then
        // decompressed data is in tmpDecompressBuffer.data(). However, these
        // are char pointers, you will need to convert them into right types as
        // follows:
        if (type == DataType::Double || type == DataType::DoubleComplex)
        {
            // TODO: use reinterpret_cast<double*>(dataIn) and
            // reinterpret_cast<double*>(tmpDecompressBuffer.data())
            // to read original data and decompressed data
        }
        else if (type == DataType::Float || type == DataType::FloatComplex)
        {
            // TODO: use reinterpret_cast<float*>(dataIn) and
            // reinterpret_cast<double*>(tmpDecompressBuffer.data())
            // to read original data and decompressed data
        }
        // TODO: after your algorithm is done, put the result into
        // *reinterpret_cast<double*>(bufferOut+bufferOutOffset) for your first
        // double number *reinterpret_cast<double*>(bufferOut+bufferOutOffset+8)
        // for your second double number and so on
    }

    bufferOutOffset += mgardBufferSize;
    return bufferOutOffset;
}

size_t CompressMGARDPlus::GetHeaderSize() const { return headerSize; }

size_t CompressMGARDPlus::DecompressV1(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    // Do NOT remove even if the buffer version is updated. Data might be still
    // in lagacy formats. This function must be kept for backward compatibility.
    // If a newer buffer format is implemented, create another function, e.g.
    // DecompressV2 and keep this function for decompressing lagacy data.

    // TODO: read your results here from
    // *reinterpret_cast<double*>(bufferIn) for your first double number
    // *reinterpret_cast<double*>(bufferIn+8) for your second double number and
    // so on

    size_t bufferInOffset = 32; // this number needs to be the same as the
                                // memory space you reserved in Operate()

    CompressMGARD mgard(m_Parameters);
    size_t sizeOut =
        mgard.InverseOperate(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOut);

    // TODO: the regular decompressed buffer is in dataOut, with the size of
    // sizeOut. Here you may want to do your magic to change the decompressed
    // data somehow to improve its accuracy :)

    headerSize += (bufferInOffset + mgard.GetHeaderSize());
    return sizeOut;
}

size_t CompressMGARDPlus::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    size_t bufferInOffset = 1; // skip operator type
    const uint8_t bufferVersion = GetParameter<uint8_t>(bufferIn, bufferInOffset);
    bufferInOffset += 2; // skip two reserved bytes
    headerSize = bufferInOffset;

    if (bufferVersion == 1)
    {
        return DecompressV1(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOut);
    }
    else if (bufferVersion == 2)
    {
        // TODO: if a Version 2 mgard buffer is being implemented, put it here
        // and keep the DecompressV1 routine for backward compatibility
    }
    else
    {
        helper::Throw<std::runtime_error>("Operator", "CompressMGARDPlus", "InverseOperate",
                                          "invalid mgard buffer version");
    }

    return 0;
}

bool CompressMGARDPlus::IsDataTypeValid(const DataType type) const
{
    if (type == DataType::Double || type == DataType::Float || type == DataType::DoubleComplex ||
        type == DataType::FloatComplex)
    {
        return true;
    }
    return false;
}

} // end namespace compress
} // end namespace core
} // end namespace adios2
