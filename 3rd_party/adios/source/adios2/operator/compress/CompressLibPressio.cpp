/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressLibPressio.cpp
 *
 *  Created on: Tue Apr 13, 2021
 *      Author: Robert Underwood robertu@g.clemson.edu
 */

#include "CompressLibPressio.h"

#include <cmath>     //std::ceil
#include <ios>       //std::ios_base::failure
#include <stdexcept> //std::invalid_argument

extern "C" {
#include "libpressio.h"
}

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace core
{
namespace compress
{

static pressio_dtype adios_to_libpressio_dtype(DataType var_type)
{
    if (var_type == helper::GetDataType<float>())
    {
        return pressio_float_dtype;
    }
    if (var_type == helper::GetDataType<double>())
    {
        return pressio_double_dtype;
    }
    if (var_type == helper::GetDataType<int8_t>())
    {
        return pressio_int8_dtype;
    }
    if (var_type == helper::GetDataType<int16_t>())
    {
        return pressio_int16_dtype;
    }
    if (var_type == helper::GetDataType<int32_t>())
    {
        return pressio_int32_dtype;
    }
    if (var_type == helper::GetDataType<int64_t>())
    {
        return pressio_int64_dtype;
    }
    if (var_type == helper::GetDataType<uint8_t>())
    {
        return pressio_uint8_dtype;
    }
    if (var_type == helper::GetDataType<uint16_t>())
    {
        return pressio_uint16_dtype;
    }
    if (var_type == helper::GetDataType<uint32_t>())
    {
        return pressio_uint32_dtype;
    }
    if (var_type == helper::GetDataType<uint64_t>())
    {
        return pressio_uint64_dtype;
    }
    helper::Throw<std::runtime_error>("Operator", "CompressLibPressio", "adios_to_libpressio_dtype",
                                      "unexpected datatype");
    return pressio_byte_dtype;
}

static std::vector<size_t> adios_to_libpressio_dims(Dims const &dims)
{
    std::vector<size_t> lp_dims(std::rbegin(dims), std::rend(dims));
    return lp_dims;
}

struct pressio_param
{
    enum class type
    {
        early,
        late,
        unset,
        malformed
    } type = type::unset;
    bool has_index = false;
    size_t index = 0;
    std::string name;
};

pressio_param parse_adios_config_entry(std::string const &key)
{
    pressio_param p;
    static const std::string early_config_prefix = "early:";
    static const std::string compressor_config_prefix = "config:";
    size_t current = 0;

    try
    {
        if (key.find(early_config_prefix) == 0)
        {
            p.type = pressio_param::type::early;
            current += early_config_prefix.size();
        }
        else if (key.find(compressor_config_prefix) == 0)
        {
            p.type = pressio_param::type::late;
            current += compressor_config_prefix.size();
        }
        else
        {
            p.type = pressio_param::type::unset;
        }
        if (current != 0)
        {
            if (key.at(current) == '[')
            {
                // we have an array entry
                p.has_index = true;
                auto digit_len = key.find_first_of(']', current);
                if (digit_len == std::string::npos)
                {
                    helper::Throw<std::invalid_argument>("Operator", "CompressLibPressio",
                                                         "parse_adios_config_entry",
                                                         "invalid substr");
                }
                p.index = stoll(key.substr(current + 1, digit_len - (current + 1)));
                current = digit_len + 2;
                if (key.at(digit_len + 1) != ':')
                {
                    helper::Throw<std::invalid_argument>("Operator", "CompressLibPressio",
                                                         "parse_adios_config_entry",
                                                         "missing expected :");
                }
            }
            else
            {
                // we have a scalar entry
                p.has_index = false;
            }

            p.name = key.substr(current);
        }
    }
    catch (std::invalid_argument &)
    {
        p.type = pressio_param::type::malformed;
    }
    catch (std::out_of_range &)
    {
        p.type = pressio_param::type::malformed;
    }
    return p;
}

static pressio_compressor *adios_to_libpressio_compressor(Params const &params)
{
    pressio *instance = pressio_instance();
    auto compressor_it = params.find("compressor_id");
    if (compressor_it != params.end())
    {
        pressio_compressor *compressor =
            pressio_get_compressor(instance, compressor_it->second.c_str());
        pressio_release(instance);
        if (compressor == nullptr)
        {
            helper::Throw<std::runtime_error>("Operator", "CompressLibPressio",
                                              "adios_to_libpressio_compressor",
                                              "compressor unavailable: " + compressor_it->second);
        }

        // adios parameters have unique names and must have string type
        // use the syntax early:foobar to set parameter foobar "early"
        // use the syntax config:foobar to set parameter foobar "late"
        // we need unique names for array entries use
        // early_config:[(\d+)]:foobar to create an array
        std::map<std::string, std::map<size_t, std::string>> early, late;
        for (auto const &param : params)
        {
            auto parsed = parse_adios_config_entry(param.first);
            std::map<std::string, std::map<size_t, std::string>> *config = nullptr;
            switch (parsed.type)
            {
            case pressio_param::type::early:
                config = &early;
                break;
            case pressio_param::type::late:
                config = &late;
                break;
            case pressio_param::type::unset:
                continue;
            case pressio_param::type::malformed:
                pressio_compressor_release(compressor);
                helper::Throw<std::runtime_error>("Operator", "CompressLibPressio",
                                                  "adios_to_libpressio_compressor",
                                                  "malformed parameter name " + param.first);
            }

            if (parsed.has_index)
            {
                (*config)[parsed.name][parsed.index] = param.second;
            }
            else
            {
                (*config)[parsed.name][static_cast<size_t>(0)] = param.second;
            }
        }

        // convert early
        pressio_options *early_opts = pressio_options_new();
        for (auto const &entry : early)
        {
            if (entry.second.size() == 1)
            {
                pressio_options_set_string(
                    early_opts, entry.first.c_str(),
                    entry.second.find(static_cast<size_t>(0))->second.c_str());
            }
            else
            {
                std::vector<const char *> entries;
                std::transform(entry.second.cbegin(), entry.second.cend(),
                               std::back_inserter(entries),
                               [](decltype(*entry.second.begin()) s) { return s.second.c_str(); });
                pressio_options_set_strings(early_opts, entry.first.c_str(), entries.size(),
                                            entries.data());
            }
        }

        pressio_compressor_set_options(compressor, early_opts);
        pressio_options_free(early_opts);

        pressio_options *compressor_options = pressio_compressor_get_options(compressor);

        for (auto const &entry : late)
        {
            pressio_option *option = nullptr;
            if (entry.second.size() == 1)
            {
                option = pressio_option_new_string(
                    entry.second.find(static_cast<size_t>(0))->second.c_str());
            }
            else
            {
                std::vector<const char *> entries;
                std::transform(entry.second.begin(), entry.second.end(),
                               std::back_inserter(entries),
                               [](decltype(*entry.second.begin()) s) { return s.second.c_str(); });
                option = pressio_option_new_strings(entries.data(), entries.size());
            }

            switch (pressio_options_cast_set(compressor_options, entry.first.c_str(), option,
                                             pressio_conversion_special))
            {
            case pressio_options_key_set:
                break;
            case pressio_options_key_exists:
                pressio_options_free(compressor_options);
                pressio_compressor_release(compressor);
                helper::Throw<std::runtime_error>("Operator", "CompressLibPressio",
                                                  "adios_to_libpressio_compressor",
                                                  "enable to convert " + entry.first);
            case pressio_options_key_does_not_exist:
                pressio_options_free(compressor_options);
                pressio_compressor_release(compressor);
                helper::Throw<std::runtime_error>("Operator", "CompressLibPressio",
                                                  "adios_to_libpressio_compressor",
                                                  "unexpected option " + entry.first);
            }
            pressio_option_free(option);
        }
        pressio_compressor_set_options(compressor, compressor_options);
        pressio_options_free(compressor_options);

        return compressor;
    }
    helper::Throw<std::runtime_error>("Operator", "CompressLibPressio",
                                      "adios_to_libpressio_compressor",
                                      "missing required \"compressor_id\" setting");
    return 0;
}

CompressLibPressio::CompressLibPressio(const Params &parameters)
: Operator("libpressio", COMPRESS_LIBPRESSIO, "compress", parameters)
{
}
size_t CompressLibPressio::Operate(const char *dataIn, const Dims &blockStart,
                                   const Dims &blockCount, const DataType type, char *bufferOut)
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
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(pressio_major_version()));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(pressio_minor_version()));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(pressio_patch_version()));
    PutParameters(bufferOut, bufferOutOffset, m_Parameters);
    // zfp V1 metadata end

    auto inputs_dims = adios_to_libpressio_dims(blockCount);
    pressio_data *input_buf =
        pressio_data_new_nonowning(adios_to_libpressio_dtype(type), const_cast<char *>(dataIn),
                                   inputs_dims.size(), inputs_dims.data());
    pressio_data *output_buf = pressio_data_new_empty(pressio_byte_dtype, 0, nullptr);
    pressio_compressor *compressor = nullptr;
    try
    {
        compressor = adios_to_libpressio_compressor(m_Parameters);
    }
    catch (std::exception &e)
    {
        pressio_data_free(input_buf);
        pressio_data_free(output_buf);
        helper::Throw<std::runtime_error>("Operator", "CompressLibPressio", "Operate",
                                          "adios_to_libpressio_compressor failed");
    }

    if (pressio_compressor_compress(compressor, input_buf, output_buf) != 0)
    {
        pressio_data_free(input_buf);
        pressio_data_free(output_buf);
        helper::Throw<std::runtime_error>(
            "Operator", "CompressLibPressio", "Operate",
            "pressio_compressor_compress: " +
                std::string(pressio_compressor_error_msg(compressor)));
    }

    size_t size_in_bytes = 0;
    void *bytes = pressio_data_ptr(output_buf, &size_in_bytes);
    memcpy(bufferOut + bufferOutOffset, bytes, size_in_bytes);
    bufferOutOffset += size_in_bytes;

    pressio_data_free(input_buf);
    pressio_data_free(output_buf);

    return bufferOutOffset;
}

size_t CompressLibPressio::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
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
        // TODO: if a Version 2 LibPressio buffer is being implemented, put it
        // here and keep the DecompressV1 routine for backward compatibility
    }
    else
    {
        helper::Throw<std::runtime_error>("Operator", "CompressLibPressio", "InverseOperate",
                                          "invalid LibPressio buffer version");
    }

    return 0;
}

bool CompressLibPressio::IsDataTypeValid(const DataType type) const
{
    if (type == DataType::Int8 || type == DataType::UInt8 || type == DataType::Int16 ||
        type == DataType::UInt16 || type == DataType::Int32 || type == DataType::UInt32 ||
        type == DataType::Int64 || type == DataType::UInt64 || type == DataType::Float ||
        type == DataType::Double)
    {
        return true;
    }
    return false;
}

size_t CompressLibPressio::DecompressV1(const char *bufferIn, const size_t sizeIn, char *dataOut)
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
    m_VersionInfo = " Data is compressed using LibPressio Version " +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) +
                    ". Please make sure a compatible version is used for decompression.";
    const Params parameters = GetParameters(bufferIn, bufferInOffset);

    std::vector<size_t> dims = adios_to_libpressio_dims(blockCount);
    pressio_data *output_buf =
        pressio_data_new_owning(adios_to_libpressio_dtype(type), dims.size(), dims.data());

    size_t newSizeIn = sizeIn - bufferInOffset;
    pressio_data *input_buf = pressio_data_new_nonowning(
        pressio_byte_dtype, const_cast<char *>(bufferIn + bufferInOffset), 1, &newSizeIn);

    pressio_compressor *compressor = nullptr;
    try
    {
        compressor = adios_to_libpressio_compressor(parameters);
    }
    catch (std::exception &)
    {
        pressio_data_free(input_buf);
        pressio_data_free(output_buf);
        helper::Throw<std::runtime_error>("Operator", "CompressLibPressio", "DecompressV1",
                                          m_VersionInfo);
    }

    if (pressio_compressor_decompress(compressor, input_buf, output_buf) != 0)
    {
        pressio_data_free(input_buf);
        pressio_data_free(output_buf);
        helper::Throw<std::runtime_error>("Operator", "CompressLibPressio", "DecompressV1",
                                          std::string("pressio_compressor_decompress: ") +
                                              pressio_compressor_error_msg(compressor) +
                                              m_VersionInfo);
    }

    size_t size_in_bytes = 0;
    void *output = pressio_data_ptr(output_buf, &size_in_bytes);
    std::memcpy(dataOut, output, size_in_bytes);

    pressio_data_free(input_buf);
    pressio_data_free(output_buf);
    return size_in_bytes;
}

} // end namespace compress
} // end namespace core
} // end namespace adios2
