/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressZFP.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include "CompressZFP.h"
#include "adios2/helper/adiosFunctions.h"
#include <sstream>
#include <zfp.h>

/* CMake will make sure zfp >= 0.5.3
   ZFP will default to SERIAL if CUDA is not available */
#ifdef ADIOS2_HAVE_ZFP_CUDA
#define ZFP_DEFAULT_EXECUTION_POLICY zfp_exec_cuda
#else
#define ZFP_DEFAULT_EXECUTION_POLICY zfp_exec_serial
#endif

namespace adios2
{
namespace core
{
namespace compress
{

/**
 * Constructor Zfp zfp_field based on input information around the data
 * pointer
 * @param data
 * @param shape
 * @param type
 * @return zfp_field*
 */
zfp_field *GetZFPField(const char *data, const Dims &shape, DataType type);

/**
 * Returns Zfp supported zfp_type based on adios string type
 * @param type adios type as string, see GetDataType<T> in
 * helper/adiosType.inl
 * @return zfp_type
 */
zfp_type GetZfpType(DataType type);

zfp_stream *GetZFPStream(const Dims &dimensions, DataType type, const Params &parameters);

CompressZFP::CompressZFP(const Params &parameters)
: Operator("zfp", COMPRESS_ZFP, "compress", parameters)
{
}

size_t CompressZFP::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                            const DataType type, char *bufferOut)
{
    const uint8_t bufferVersion = 1;
    size_t bufferOutOffset = 0;

    MakeCommonHeader(bufferOut, bufferOutOffset, bufferVersion);

    const size_t ndims = blockCount.size();

    // zfp V1 metadata
    PutParameter(bufferOut, bufferOutOffset, ndims);
    for (const auto &d : blockCount)
    {
        PutParameter(bufferOut, bufferOutOffset, d);
    }
    PutParameter(bufferOut, bufferOutOffset, type);
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(ZFP_VERSION_MAJOR));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(ZFP_VERSION_MINOR));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(ZFP_VERSION_PATCH));
    PutParameters(bufferOut, bufferOutOffset, m_Parameters);
    // zfp V1 metadata end

    Dims convertedDims = ConvertDims(blockCount, type, 3);

    zfp_field *field = GetZFPField(dataIn, convertedDims, type);
    zfp_stream *stream = GetZFPStream(convertedDims, type, m_Parameters);

    size_t maxSize = zfp_stream_maximum_size(stream, field);
    // associate bitstream
    bitstream *bitstream = stream_open(bufferOut + bufferOutOffset, maxSize);
    zfp_stream_set_bit_stream(stream, bitstream);
    zfp_stream_rewind(stream);

    size_t sizeOut = zfp_compress(stream, field);

    if (sizeOut == 0)
    {
        helper::Throw<std::runtime_error>("Operator", "CompressZFP", "Operate(Compress)",
                                          "zfp failed, compressed buffer size is 0");
    }

    bufferOutOffset += sizeOut;

    zfp_field_free(field);
    zfp_stream_close(stream);
    stream_close(bitstream);

    return bufferOutOffset;
}

size_t CompressZFP::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    size_t bufferInOffset = 1; // skip operator type
    const uint8_t bufferVersion = GetParameter<uint8_t>(bufferIn, bufferInOffset);
    bufferInOffset += 2; // skip two reserved bytes

    if (bufferVersion == 1)
    {
        return DecompressV1(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOut);
    }
    else if (bufferVersion == 2)
    {
        // TODO: if a Version 2 zfp buffer is being implemented, put it here
        // and keep the DecompressV1 routine for backward compatibility
    }
    else
    {
        helper::Throw<std::runtime_error>("Operator", "CompressZFP", "InverseOperate",
                                          "invalid zfp buffer version" +
                                              std::to_string(bufferVersion));
    }

    return 0;
}

bool CompressZFP::IsDataTypeValid(const DataType type) const
{
    if (type == DataType::Float || type == DataType::Double || type == DataType::FloatComplex ||
        type == DataType::DoubleComplex || type == DataType::Int32 || type == DataType::Int64)
    {
        return true;
    }
    return false;
}

// PRIVATE

size_t CompressZFP::DecompressV1(const char *bufferIn, const size_t sizeIn, char *dataOut)
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
    this->m_VersionInfo = " Data is compressed using ZFP Version " +
                          std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                          std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                          std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) +
                          ". Please make sure a compatible version is used for decompression.";
    const Params parameters = GetParameters(bufferIn, bufferInOffset);

    Dims convertedDims = ConvertDims(blockCount, type, 3);

    zfp_field *field = nullptr;
    zfp_stream *stream = nullptr;

    field = GetZFPField(dataOut, convertedDims, type);
    stream = GetZFPStream(convertedDims, type, parameters);

    // associate bitstream
    bitstream *bitstream =
        stream_open(const_cast<char *>(bufferIn + bufferInOffset), sizeIn - bufferInOffset);
    zfp_stream_set_bit_stream(stream, bitstream);
    zfp_stream_rewind(stream);

    int status = zfp_decompress(stream, field);

    if (!status)
    {
        helper::Throw<std::runtime_error>("Operator", "CompressZFP", "DecompressV1",
                                          "zfp failed with status " + std::to_string(status));
    }

    zfp_field_free(field);
    zfp_stream_close(stream);
    stream_close(bitstream);

    return helper::GetTotalSize(convertedDims, helper::GetDataTypeSize(type));
}

zfp_type GetZfpType(DataType type)
{
    zfp_type zfpType = zfp_type_none;

    if (type == helper::GetDataType<double>())
    {
        zfpType = zfp_type_double;
    }
    else if (type == helper::GetDataType<float>())
    {
        zfpType = zfp_type_float;
    }
    else if (type == helper::GetDataType<int64_t>())
    {
        zfpType = zfp_type_int64;
    }
    else if (type == helper::GetDataType<int32_t>())
    {
        zfpType = zfp_type_int32;
    }
    else if (type == helper::GetDataType<std::complex<float>>())
    {
        zfpType = zfp_type_float;
    }
    else if (type == helper::GetDataType<std::complex<double>>())
    {
        zfpType = zfp_type_double;
    }
    else
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressZFP", "GetZfpType",
                                             "invalid data type " + ToString(type));
    }

    return zfpType;
}

// Free static functions

zfp_field *GetZFPField(const char *data, const Dims &dimensions, DataType type)
{
    zfp_type zfpType = GetZfpType(type);
    zfp_field *field = nullptr;

    if (dimensions.size() == 1)
    {
        field = zfp_field_1d(const_cast<char *>(data), zfpType, dimensions[0]);
    }
    else if (dimensions.size() == 2)
    {
        field = zfp_field_2d(const_cast<char *>(data), zfpType, dimensions[0], dimensions[1]);
    }
    else if (dimensions.size() == 3)
    {
        field = zfp_field_3d(const_cast<char *>(data), zfpType, dimensions[0], dimensions[1],
                             dimensions[2]);
    }
    else
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressZFP", "GetZfpField",
                                             "zfp does not support " +
                                                 std::to_string(dimensions.size()) + "D data");
    }

    if (field == nullptr)
    {
        helper::Throw<std::runtime_error>("Operator", "CompressZFP", "GetZfpField",
                                          "zfp failed to make field for" +
                                              std::to_string(dimensions.size()) + "D data in " +
                                              ToString(type));
    }

    return field;
}

zfp_stream *GetZFPStream(const Dims &dimensions, DataType type, const Params &parameters)
{
    zfp_stream *stream = zfp_stream_open(NULL);
    bool isSerial = true;

    auto itAccuracy = parameters.find("accuracy");
    const bool hasAccuracy = itAccuracy != parameters.end();

    auto itRate = parameters.find("rate");
    const bool hasRate = itRate != parameters.end();

    auto itPrecision = parameters.find("precision");
    const bool hasPrecision = itPrecision != parameters.end();

    auto itBackend = parameters.find("backend");
    const bool hasBackend = itBackend != parameters.end();

    zfp_stream_set_execution(stream, ZFP_DEFAULT_EXECUTION_POLICY);
    isSerial = ZFP_DEFAULT_EXECUTION_POLICY == zfp_exec_serial;

#ifdef ADIOS2_HAVE_ZFP_CUDA
    if (hasAccuracy)
    {
        helper::Throw<std::runtime_error>(
            "Operator", "CompressZFP", "GetZfpField",
            "The CUDA backend in ZFP cannot use the 'accuracy' parameter");
    }
    if (hasPrecision)
    {
        helper::Throw<std::runtime_error>(
            "Operator", "CompressZFP", "GetZfpField",
            "The CUDA backend in ZFP cannot use the 'precisiom' parameter");
    }
#endif

    if (hasBackend)
    {
        auto policy = ZFP_DEFAULT_EXECUTION_POLICY;
        const auto backend = itBackend->second;

        if (backend == "serial")
        {
            policy = zfp_exec_serial;
            isSerial = true;
        }
        else if (backend == "omp")
        {
            policy = zfp_exec_omp;
        }
#ifdef ADIOS2_HAVE_ZFP_CUDA
        else if (backend == "cuda")
        {
            policy = zfp_exec_cuda;
        }
#endif

        zfp_stream_set_execution(stream, policy);
    }

    if ((hasAccuracy && hasPrecision) || (hasAccuracy && hasRate) || (hasPrecision && hasRate) ||
        (!hasAccuracy && !hasRate && !hasPrecision && !isSerial))
    {
        std::ostringstream oss;
        oss << std::endl
            << "ERROR: ZFP:"
               " 'accuracy'|'rate'|'precision' params are mutualy exclusive: ";

        for (auto &p : parameters)
        {
            oss << "(" << p.first << ", " << p.second << ").";
        }
        helper::Throw<std::runtime_error>("Operator", "CompressZFP", "GetZfpField", oss.str());
    }

    if (hasAccuracy)
    {
        const double accuracy = helper::StringTo<double>(
            itAccuracy->second, "setting 'accuracy' in call to CompressZfp\n");

        zfp_stream_set_accuracy(stream, accuracy);
    }
    else if (hasRate)
    {
        const double rate =
            helper::StringTo<double>(itRate->second, "setting 'rate' in call to CompressZfp\n");
        // TODO support last argument write random access?
        zfp_stream_set_rate(stream, rate, GetZfpType(type),
                            static_cast<unsigned int>(dimensions.size()), 0);
    }
    else if (hasPrecision)
    {
        const unsigned int precision = static_cast<unsigned int>(helper::StringTo<uint32_t>(
            itPrecision->second, "setting 'precision' in call to CompressZfp\n"));
        zfp_stream_set_precision(stream, precision);
    }

    return stream;
}

} // end namespace compress
} // end namespace core
} // end namespace adios2
