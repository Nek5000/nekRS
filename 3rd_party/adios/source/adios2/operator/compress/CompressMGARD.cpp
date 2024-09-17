/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressMGARD.cpp :
 *
 *  Created on: Dec 1, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "CompressMGARD.h"
#include "CompressNull.h"
#include "adios2/helper/adiosFunctions.h"
#include <cstring>
#include <mgard/MGARDConfig.hpp>
#include <mgard/compress_x.hpp>

namespace adios2
{
namespace core
{
namespace compress
{

CompressMGARD::CompressMGARD(const Params &parameters)
: Operator("mgard", COMPRESS_MGARD, "compress", parameters)
{
}

size_t CompressMGARD::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                              const DataType type, char *bufferOut)
{
    const uint8_t bufferVersion = 1;
    size_t bufferOutOffset = 0;

    MakeCommonHeader(bufferOut, bufferOutOffset, bufferVersion);

    Dims convertedDims = ConvertDims(blockCount, type, 3);

    const size_t ndims = convertedDims.size();
    if (ndims > 5)
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressMGARD", "Operate",
                                             "MGARD does not support data in " +
                                                 std::to_string(ndims) + " dimensions");
    }

    // mgard V1 metadata
    PutParameter(bufferOut, bufferOutOffset, ndims);
    for (const auto &d : convertedDims)
    {
        PutParameter(bufferOut, bufferOutOffset, d);
    }
    PutParameter(bufferOut, bufferOutOffset, type);
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(MGARD_VERSION_MAJOR));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(MGARD_VERSION_MINOR));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(MGARD_VERSION_PATCH));
    // mgard V1 metadata end

    // set type
    mgard_x::data_type mgardType;
    if (type == helper::GetDataType<float>())
    {
        mgardType = mgard_x::data_type::Float;
    }
    else if (type == helper::GetDataType<double>())
    {
        mgardType = mgard_x::data_type::Double;
    }
    else if (type == helper::GetDataType<std::complex<float>>())
    {
        mgardType = mgard_x::data_type::Float;
    }
    else if (type == helper::GetDataType<std::complex<double>>())
    {
        mgardType = mgard_x::data_type::Double;
    }
    else
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressMGARD", "Operate",
                                             "MGARD only supports float and double types");
    }
    // set type end

    // set mgard style dim info
    mgard_x::DIM mgardDim = ndims;
    std::vector<mgard_x::SIZE> mgardCount;
    for (const auto &c : convertedDims)
    {
        mgardCount.push_back(c);
    }
    // set mgard style dim info end

    // Parameters
    bool hasTolerance = false;
    double tolerance = 0.0;
    double s = 0.0;
    auto errorBoundType = mgard_x::error_bound_type::REL;

    // input size under this bound will not compress
    size_t thresholdSize = 100000;

    auto itThreshold = m_Parameters.find("threshold");
    if (itThreshold != m_Parameters.end())
    {
        thresholdSize = std::stod(itThreshold->second);
    }
    auto itAccuracy = m_Parameters.find("accuracy");
    if (itAccuracy != m_Parameters.end())
    {
        tolerance = std::stod(itAccuracy->second);
        hasTolerance = true;
    }
    auto itTolerance = m_Parameters.find("tolerance");
    if (itTolerance != m_Parameters.end())
    {
        tolerance = std::stod(itTolerance->second);
        hasTolerance = true;
    }
    if (!hasTolerance)
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressMGARD", "Operate",
                                             "missing mandatory parameter tolerance / accuracy");
    }
    auto itSParameter = m_Parameters.find("s");
    if (itSParameter != m_Parameters.end())
    {
        s = std::stod(itSParameter->second);
    }
    auto itMode = m_Parameters.find("mode");
    if (itMode != m_Parameters.end())
    {
        if (itMode->second == "ABS")
        {
            errorBoundType = mgard_x::error_bound_type::ABS;
        }
        else if (itMode->second == "REL")
        {
            errorBoundType = mgard_x::error_bound_type::REL;
        }
    }

    // let mgard know the output buffer size
    size_t sizeOut = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(type));

    if (sizeOut < thresholdSize)
    {
        /* disable compression and add marker in the header*/
        PutParameter(bufferOut, bufferOutOffset, false);
        headerSize = bufferOutOffset;
        return 0;
    }

    PutParameter(bufferOut, bufferOutOffset, true);
    void *compressedData = bufferOut + bufferOutOffset;
    mgard_x::compress(mgardDim, mgardType, mgardCount, tolerance, s, errorBoundType, dataIn,
                      compressedData, sizeOut, true);
    bufferOutOffset += sizeOut;

    return bufferOutOffset;
}

size_t CompressMGARD::GetHeaderSize() const { return headerSize; }

size_t CompressMGARD::DecompressV1(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    // Do NOT remove even if the buffer version is updated. Data might be still
    // in lagacy formats. This function must be kept for backward compatibility.
    // If a newer buffer format is implemented, create another function, e.g.
    // DecompressV2 and keep this function for decompressing lagacy data.

    size_t bufferInOffset = 0;

    const size_t ndims = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    Dims blockCount(ndims);
    for (size_t i = 0; i < ndims; ++i)
    {
        blockCount[i] = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    }
    const DataType type = GetParameter<DataType>(bufferIn, bufferInOffset);
    m_VersionInfo = " Data is compressed using MGARD Version " +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) +
                    ". Please make sure a compatible version is used for decompression.";

    const bool isCompressed = GetParameter<bool>(bufferIn, bufferInOffset);

    size_t sizeOut = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(type));

    if (type == DataType::FloatComplex || type == DataType::DoubleComplex)
    {
        sizeOut /= 2;
    }

    if (isCompressed)
    {
        try
        {
            void *dataOutVoid = dataOut;
            mgard_x::decompress(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOutVoid,
                                true);
        }
        catch (...)
        {
            helper::Throw<std::runtime_error>("Operator", "CompressMGARD", "DecompressV1",
                                              m_VersionInfo);
        }
        return sizeOut;
    }

    headerSize += bufferInOffset;
    return 0;
}

size_t CompressMGARD::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
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
        helper::Throw<std::runtime_error>("Operator", "CompressMGARD", "InverseOperate",
                                          "invalid mgard buffer version");
    }

    return 0;
}

bool CompressMGARD::IsDataTypeValid(const DataType type) const
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
