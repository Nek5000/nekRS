/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * it->first ==ompanying file Copyright.txt for details.
 *
 * CompressSZ.cpp
 *
 *  Created on: Jul 24, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "CompressSZ.h"
#include "CompressNull.h"
#include "adios2/helper/adiosFunctions.h"
#include <cmath>     //std::ceil
#include <ios>       //std::ios_base::failure
#include <stdexcept> //std::invalid_argument
extern "C" {
#include <sz.h>
}

namespace adios2
{
namespace core
{
namespace compress
{

std::mutex CompressSZ::m_Mutex;

CompressSZ::CompressSZ(const Params &parameters)
: Operator("sz", COMPRESS_SZ, "compress", parameters)
{
}

size_t CompressSZ::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                           const DataType varType, char *bufferOut)
{
    const uint8_t bufferVersion = 2;
    size_t bufferOutOffset = 0;

    MakeCommonHeader(bufferOut, bufferOutOffset, bufferVersion);

    Dims convertedDims = ConvertDims(blockCount, varType, 5);
    const size_t ndims = convertedDims.size();

    // sz V2 metadata
    PutParameter(bufferOut, bufferOutOffset, ndims);
    for (const auto &d : convertedDims)
    {
        PutParameter(bufferOut, bufferOutOffset, d);
    }
    PutParameter(bufferOut, bufferOutOffset, varType);
    for (uint8_t i = 0; i < 4; ++i)
    {
        PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(versionNumber[i]));
    }
    // sz V2 metadata end

    sz_params sz;
    memset(&sz, 0, sizeof(sz_params));
    sz.max_quant_intervals = 65536;
    sz.quantization_intervals = 0;
    //    sz.dataEndianType = LITTLE_ENDIAN_DATA;
    //    sz.sysEndianType = LITTLE_ENDIAN_DATA;
    sz.sol_ID = SZ;
    // sz.layers = 1;
    sz.sampleDistance = 100;
    sz.predThreshold = 0.99;
    //    sz.offset = 0;
    sz.szMode = SZ_BEST_COMPRESSION; // SZ_BEST_SPEED; //SZ_BEST_COMPRESSION;
    sz.gzipMode = 1;
    sz.errorBoundMode = ABS;
    sz.absErrBound = 1E-4;
    sz.relBoundRatio = 1E-3;
    sz.psnr = 80.0;
    sz.pw_relBoundRatio = 1E-5;
    sz.segment_size = static_cast<int>(std::pow(5., static_cast<double>(ndims)));
    sz.pwr_type = SZ_PWR_MIN_TYPE;

    /* SZ parameters */
    int use_configfile = 0;
    std::string sz_configfile = "sz.config";

    Params::const_iterator it;
    for (it = m_Parameters.begin(); it != m_Parameters.end(); it++)
    {
        if (it->first == "init")
        {
            use_configfile = 1;
            sz_configfile = std::string(it->second);
        }
        else if (it->first == "max_quant_intervals")
        {
            sz.max_quant_intervals = std::stoi(it->second);
        }
        else if (it->first == "quantization_intervals")
        {
            sz.quantization_intervals = std::stoi(it->second);
        }
        else if (it->first == "sol_ID")
        {
            sz.sol_ID = std::stoi(it->second);
        }
        else if (it->first == "sampleDistance")
        {
            sz.sampleDistance = std::stoi(it->second);
        }
        else if (it->first == "predThreshold")
        {
            sz.predThreshold = std::stof(it->second);
        }
        else if (it->first == "szMode")
        {
            int szMode = SZ_BEST_SPEED;
            if (it->second == "SZ_BEST_SPEED")
            {
                szMode = SZ_BEST_SPEED;
            }
            else if (it->second == "SZ_BEST_COMPRESSION")
            {
                szMode = SZ_BEST_COMPRESSION;
            }
            else if (it->second == "SZ_DEFAULT_COMPRESSION")
            {
                szMode = SZ_DEFAULT_COMPRESSION;
            }
            else
            {
                helper::Throw<std::invalid_argument>(
                    "Operator", "CompressSZ", "Operate",
                    "Parameter szMode must be SZ_BEST_SPEED, "
                    "SZ_BEST_COMPRESSION or SZ_DEFAULT_COMPRESSION");
            }
            sz.szMode = szMode;
        }
        else if (it->first == "gzipMode")
        {
            sz.gzipMode = std::stoi(it->second);
        }
        else if (it->first == "errorBoundMode")
        {
            int errorBoundMode = ABS;
            if (it->second == "ABS")
            {
                errorBoundMode = ABS;
            }
            else if (it->second == "REL")
            {
                errorBoundMode = REL;
            }
            else if (it->second == "ABS_AND_REL")
            {
                errorBoundMode = ABS_AND_REL;
            }
            else if (it->second == "ABS_OR_REL")
            {
                errorBoundMode = ABS_OR_REL;
            }
            else if (it->second == "PW_REL")
            {
                errorBoundMode = PW_REL;
            }
            else
            {
                helper::Throw<std::invalid_argument>(
                    "Operator", "CompressSZ", "Operate",
                    "Parameter errorBoundMode must be ABS, REL, ABS_AND_REL, "
                    "ABS_OR_REL or PW_REL");
            }
            sz.errorBoundMode = errorBoundMode;
        }
        else if (it->first == "absErrBound")
        {
            sz.absErrBound = std::stof(it->second);
        }
        else if (it->first == "relBoundRatio")
        {
            sz.relBoundRatio = std::stof(it->second);
        }
        else if (it->first == "pw_relBoundRatio")
        {
            sz.pw_relBoundRatio = std::stof(it->second);
        }
        else if (it->first == "segment_size")
        {
            sz.segment_size = std::stoi(it->second);
        }
        else if (it->first == "pwr_type")
        {
            int pwr_type = SZ_PWR_MIN_TYPE;
            if ((it->second == "MIN") || (it->second == "SZ_PWR_MIN_TYPE"))
            {
                pwr_type = SZ_PWR_MIN_TYPE;
            }
            else if ((it->second == "AVG") || (it->second == "SZ_PWR_AVG_TYPE"))
            {
                pwr_type = SZ_PWR_AVG_TYPE;
            }
            else if ((it->second == "MAX") || (it->second == "SZ_PWR_MAX_TYPE"))
            {
                pwr_type = SZ_PWR_MAX_TYPE;
            }
            else
            {
                helper::Throw<std::invalid_argument>(
                    "Operator", "CompressSZ", "Operate",
                    "Parameter pwr_type must be MIN, SZ_PWR_MIN_TYPE, AVG, "
                    "SZ_PWR_AVG_TYPE, MAX or SZ_PWR_MAX_TYPE");
            }
            sz.pwr_type = pwr_type;
        }
        else if ((it->first == "abs") || (it->first == "absolute") || (it->first == "accuracy"))
        {
            sz.errorBoundMode = ABS;
            sz.absErrBound = std::stod(it->second);
        }
        else if ((it->first == "rel") || (it->first == "relative"))
        {
            sz.errorBoundMode = REL;
            sz.relBoundRatio = std::stof(it->second);
        }
        else if ((it->first == "pw") || (it->first == "pwr") || (it->first == "pwrel") ||
                 (it->first == "pwrelative"))
        {
            sz.errorBoundMode = PW_REL;
            sz.pw_relBoundRatio = std::stof(it->second);
        }
        else if ((it->first == "zchecker") || (it->first == "zcheck") ||
                 (it->first == "z-checker") || (it->first == "z-check"))
        {
            // TODO:
            // Z-checker is not currently implemented
            // use_zchecker = (it->second == "") ? 1 : std::stof(it->second);
        }
        else
        {
            // TODO: ignoring unknown option due to Fortran bindings passing
            // empty parameter
        }
    }

    // Get type info
    int dtype = -1;
    if (varType == helper::GetDataType<double>() ||
        varType == helper::GetDataType<std::complex<double>>())
    {
        dtype = SZ_DOUBLE;
    }
    else if (varType == helper::GetDataType<float>() ||
             varType == helper::GetDataType<std::complex<float>>())
    {
        dtype = SZ_FLOAT;
    }
    else
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressSZ", "Operate",
                                             "SZ compressor only support float or double types");
    }

    convertedDims = ConvertDims(blockCount, varType, 5, true, 0);

    size_t szBufferSize;

    m_Mutex.lock();
    if (use_configfile)
    {
        SZ_Init(sz_configfile.c_str());
    }
    else
    {
        SZ_Init_Params(&sz);
    }
    auto *szBuffer =
        SZ_compress(dtype, const_cast<char *>(dataIn), &szBufferSize, convertedDims[0],
                    convertedDims[1], convertedDims[2], convertedDims[3], convertedDims[4]);
    SZ_Finalize();
    m_Mutex.unlock();

    if (bufferOutOffset + szBufferSize >
        helper::GetTotalSize(blockCount, helper::GetDataTypeSize(varType)))
    {
        CompressNull c({});
        bufferOutOffset = c.Operate(dataIn, blockStart, blockCount, varType, bufferOut);
    }
    else
    {
        std::memcpy(bufferOut + bufferOutOffset, szBuffer, szBufferSize);
        bufferOutOffset += szBufferSize;
    }
    if (szBuffer)
    {
        free(szBuffer);
    }
    return bufferOutOffset;
}

size_t CompressSZ::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    size_t bufferInOffset = 1; // skip operator type
    const uint8_t bufferVersion = GetParameter<uint8_t>(bufferIn, bufferInOffset);
    bufferInOffset += 2; // skip two reserved bytes

    if (bufferVersion == 2)
    {
        return DecompressV2(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOut);
    }
    else if (bufferVersion == 1)
    {
        return DecompressV1(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOut);
    }
    else if (bufferVersion == 3)
    {
        // TODO: if a Version 3 sz buffer is being implemented, put it here
        // and keep the DecompressV1 routine for backward compatibility
    }
    else
    {
        helper::Throw<std::runtime_error>("Operator", "CompressSZ", "InverseOperate",
                                          "invalid sz buffer version");
    }

    return 0;
}

bool CompressSZ::IsDataTypeValid(const DataType type) const
{
    if (type == DataType::Float || type == DataType::Double || type == DataType::FloatComplex ||
        type == DataType::DoubleComplex)
    {
        return true;
    }
    return false;
}

size_t CompressSZ::DecompressV1(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    // Do NOT remove even if the buffer version is updated. Data might be still
    // in lagacy formats. This function must be kept for backward compatibility.
    // If a newer buffer format is implemented, create another function, e.g.
    // DecompressV2 and keep this function for decompressing legacy data.

    size_t bufferInOffset = 0;

    const size_t ndims = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    Dims blockCount(ndims);
    for (size_t i = 0; i < ndims; ++i)
    {
        blockCount[i] = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    }
    const DataType type = GetParameter<DataType>(bufferIn, bufferInOffset);

    m_VersionInfo = " Data is compressed using SZ Version " +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) +
                    ". Please make sure a compatible version is used for decompression.";

    Dims convertedDims = ConvertDims(blockCount, type, 4, true, 1);

    // Get type info
    int dtype = 0;
    size_t dataTypeSize = 0;
    if (type == helper::GetDataType<double>() ||
        type == helper::GetDataType<std::complex<double>>())
    {
        dtype = SZ_DOUBLE;
        dataTypeSize = 8;
    }
    else if (type == helper::GetDataType<float>() ||
             type == helper::GetDataType<std::complex<float>>())
    {
        dtype = SZ_FLOAT;
        dataTypeSize = 4;
    }
    else
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressSZ", "DecompressV1",
                                             "SZ compressor only support float or double types");
    }

    const size_t dataSizeBytes = helper::GetTotalSize(convertedDims, dataTypeSize);

    m_Mutex.lock();
    void *result = SZ_decompress(
        dtype, reinterpret_cast<unsigned char *>(const_cast<char *>(bufferIn + bufferInOffset)),
        sizeIn - bufferInOffset, 0, convertedDims[0], convertedDims[1], convertedDims[2],
        convertedDims[3]);
    SZ_Finalize();
    m_Mutex.unlock();

    if (result == nullptr)
    {
        helper::Throw<std::runtime_error>("Operator", "CompressSZ", "DecompressV1", m_VersionInfo);
    }
    std::memcpy(dataOut, result, dataSizeBytes);
    free(result);
    result = nullptr;
    return dataSizeBytes;
}

size_t CompressSZ::DecompressV2(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    // Do NOT remove even if the buffer version is updated. Data might be still
    // in lagacy formats. This function must be kept for backward compatibility.
    // If a newer buffer format is implemented, create another function, e.g.
    // DecompressV3 and keep this function for decompressing legacy data.

    size_t bufferInOffset = 0;

    const size_t ndims = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    Dims blockCount(ndims);
    for (size_t i = 0; i < ndims; ++i)
    {
        blockCount[i] = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    }
    const DataType type = GetParameter<DataType>(bufferIn, bufferInOffset);

    m_VersionInfo = " Data is compressed using SZ Version " +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) +
                    ". Please make sure a compatible version is used for decompression.";

    // Get type info
    int dtype = 0;
    size_t dataTypeSize = 0;
    if (type == helper::GetDataType<double>() ||
        type == helper::GetDataType<std::complex<double>>())
    {
        dtype = SZ_DOUBLE;
        dataTypeSize = 8;
    }
    else if (type == helper::GetDataType<float>() ||
             type == helper::GetDataType<std::complex<float>>())
    {
        dtype = SZ_FLOAT;
        dataTypeSize = 4;
    }
    else
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressSZ", "DecompressV2",
                                             "SZ compressor only supports float or double types");
    }

    const size_t dataSizeBytes = helper::GetTotalSize(blockCount, dataTypeSize);

    Dims convertedDims = ConvertDims(blockCount, helper::GetDataType<float>(), 5, true, 0);

    m_Mutex.lock();
    void *result = SZ_decompress(
        dtype, reinterpret_cast<unsigned char *>(const_cast<char *>(bufferIn + bufferInOffset)),
        sizeIn - bufferInOffset, convertedDims[0], convertedDims[1], convertedDims[2],
        convertedDims[3], convertedDims[4]);
    SZ_Finalize();
    m_Mutex.unlock();

    if (result == nullptr)
    {
        helper::Throw<std::runtime_error>("Operator", "CompressSZ", "DecompressV2", m_VersionInfo);
    }
    std::memcpy(dataOut, result, dataSizeBytes);
    free(result);
    result = nullptr;
    return dataSizeBytes;
}

} // end namespace compress
} // end namespace core
} // end namespace adios2
