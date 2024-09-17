/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressPNG.cpp
 *
 *  Created on: Jun 10, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include "CompressPNG.h"

#include <cstring> // std::memset

extern "C" {
#include <png.h>
}

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace core
{
namespace compress
{

const std::map<std::string, uint32_t> CompressPNG::m_ColorTypes = {
    {"PNG_COLOR_TYPE_GRAY", PNG_COLOR_TYPE_GRAY},
    {"PNG_COLOR_TYPE_PALETTE", PNG_COLOR_TYPE_PALETTE},
    {"PNG_COLOR_TYPE_RGB", PNG_COLOR_TYPE_RGB},
    {"PNG_COLOR_TYPE_RGB_ALPHA", PNG_COLOR_TYPE_RGB_ALPHA},
    {"PNG_COLOR_TYPE_GRAY_ALPHA", PNG_COLOR_TYPE_GRAY_ALPHA},
    {"PNG_COLOR_TYPE_RGBA", PNG_COLOR_TYPE_RGBA},
    {"PNG_COLOR_TYPE_GA", PNG_COLOR_TYPE_GA}};

const std::map<std::string, std::set<uint32_t>> CompressPNG::m_BitDepths = {
    {"PNG_COLOR_TYPE_GRAY", {1, 2, 4, 8, 16}},
    {"PNG_COLOR_TYPE_PALETTE", {1, 2, 4, 8}},
    {"PNG_COLOR_TYPE_RGB", {8, 16}},
    {"PNG_COLOR_TYPE_RGB_ALPHA", {8, 16}},
    {"PNG_COLOR_TYPE_GRAY_ALPHA", {8, 16}},
    {"PNG_COLOR_TYPE_RGBA", {8, 16}},
    {"PNG_COLOR_TYPE_GA", {8, 16}}};

// PUBLIC
CompressPNG::CompressPNG(const Params &parameters)
: Operator("png", COMPRESS_PNG, "compress", parameters)
{
}

size_t CompressPNG::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                            const DataType type, char *bufferOut)
{
    size_t bufferOutOffset = 0;
    const uint8_t bufferVersion = 1;

    MakeCommonHeader(bufferOut, bufferOutOffset, bufferVersion);

    size_t paramOffset = bufferOutOffset;
    bufferOutOffset += sizeof(size_t) + 3;

    auto lf_Write = [](png_structp png_ptr, png_bytep data, png_size_t length) {
        DestInfo *pDestInfo = reinterpret_cast<DestInfo *>(png_get_io_ptr(png_ptr));
        std::memcpy(pDestInfo->BufferOut + pDestInfo->Offset, data, length);
        pDestInfo->Offset += length;
    };

    const std::size_t ndims = blockCount.size();

    if (ndims != 3 && ndims != 2)
    {
        helper::Throw<std::invalid_argument>(
            "Operator", "CompressPNG", "Operate",
            "image number of dimensions " + std::to_string(ndims) +
                " is invalid, must be 2 {height,width*bytes_per_pixel} or 3"
                " {height,width,bytes_per_pixel]} , in call to ADIOS2 PNG "
                " compression");
    }

    // defaults
    int compressionLevel = 1;
    int colorType = PNG_COLOR_TYPE_RGBA;
    int bitDepth = 8;
    std::string colorTypeStr = "PNG_COLOR_TYPE_RGBA";

    for (const auto &itParameter : m_Parameters)
    {
        const std::string key = itParameter.first;
        const std::string value = itParameter.second;

        if (key == "compression_level")
        {
            compressionLevel = static_cast<int>(
                helper::StringTo<int32_t>(value, "when setting PNG level parameter\n"));

            if (compressionLevel < 1 || compressionLevel > 9)
            {
                helper::Throw<std::invalid_argument>(
                    "Operator", "CompressPNG", "Operate",
                    "compression_level must be an "
                    "integer between 1 (less "
                    "compression, less memory) and 9 "
                    "(more compression, more memory) inclusive, in call to "
                    "ADIOS2 PNG Compress");
            }
        }
        else if (key == "color_type")
        {
            auto itColorType = m_ColorTypes.find(value);

            if (itColorType == m_ColorTypes.end())
            {
                helper::Throw<std::invalid_argument>(
                    "Operator", "CompressPNG", "Operate",
                    "invalid color_type, see PNG_COLOR_TYPE_* for "
                    "available types, in call to ADIOS2 PNG Compress");
            }

            colorTypeStr = itColorType->first;
            colorType = itColorType->second;
        }
        else if (key == "bit_depth")
        {
            bitDepth = static_cast<int>(
                helper::StringTo<int32_t>(value, "when setting PNG bit_depth parameter\n"));
        }
    }

    if (m_BitDepths.at(colorTypeStr).count(static_cast<int32_t>(bitDepth)) == 0)
    {
        helper::Throw<std::invalid_argument>(
            "Operator", "CompressPNG", "Operate",
            "bit_depth " + std::to_string(bitDepth) + " and color_type " + colorTypeStr +
                " combination is not allowed by libpng, in call to ADIOS2 PNG "
                "compression");
    }

    png_structp pngWrite =
        png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop pngInfo = png_create_info_struct(pngWrite);

    const uint32_t bytesPerPixel = ndims == 3
                                       ? static_cast<uint32_t>(blockCount[2])
                                       : static_cast<uint32_t>(helper::GetDataTypeSize(type));

    const uint32_t width = static_cast<uint32_t>(blockCount[1]);
    const uint32_t height = static_cast<uint32_t>(blockCount[0]);

    png_set_IHDR(pngWrite, pngInfo, width, height, bitDepth, colorType, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    if (setjmp(png_jmpbuf(pngWrite)))
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressPNG", "Operate",
                                             "libpng detected an error in ADIOS2 PNG Compress");
    }

    png_set_compression_level(pngWrite, compressionLevel);

    // set the rows
    std::vector<uint8_t *> rows(height);
    for (size_t r = 0; r < height; ++r)
    {
        rows[r] =
            reinterpret_cast<uint8_t *>(const_cast<char *>(dataIn)) + r * width * bytesPerPixel;
    }
    png_set_rows(pngWrite, pngInfo, rows.data());

    DestInfo destInfo;
    destInfo.BufferOut = bufferOut;
    destInfo.Offset = bufferOutOffset;

    png_set_write_fn(pngWrite, &destInfo, lf_Write, nullptr);
    png_write_png(pngWrite, pngInfo, PNG_TRANSFORM_IDENTITY, nullptr);
    png_write_end(pngWrite, pngInfo);

    // const size_t compressedSize = png_get_compression_buffer_size(pngWrite);
    png_destroy_write_struct(&pngWrite, &pngInfo);

    PutParameter(bufferOut, paramOffset, destInfo.Offset);
    PutParameter(bufferOut, paramOffset, static_cast<uint8_t>(PNG_LIBPNG_VER_MAJOR));
    PutParameter(bufferOut, paramOffset, static_cast<uint8_t>(PNG_LIBPNG_VER_MINOR));
    PutParameter(bufferOut, paramOffset, static_cast<uint8_t>(PNG_LIBPNG_VER_RELEASE));

    return destInfo.Offset;
}

size_t CompressPNG::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    size_t bufferInOffset = 1; // skip operator type
    const uint8_t bufferVersion = GetParameter<uint8_t>(bufferIn, bufferInOffset);
    bufferInOffset += 2; // skip two reserved bytes

    if (bufferVersion == 1)
    {
        // pass in the whole buffer as there is absolute positions saved in the
        // buffer to determine the offsets and lengths for batches
        return DecompressV1(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOut);
    }
    else if (bufferVersion == 2)
    {
        // TODO: if a Version 2 png buffer is being implemented, put it here
        // and keep the DecompressV1 routine for backward compatibility
    }
    else
    {
        helper::Throw<std::runtime_error>("Operator", "CompressPNG", "InverseOperate",
                                          "invalid png buffer version");
    }

    return 0;
}

bool CompressPNG::IsDataTypeValid(const DataType type) const { return true; }

size_t CompressPNG::DecompressV1(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    // Do NOT remove even if the buffer version is updated. Data might be still
    // in lagacy formats. This function must be kept for backward compatibility.
    // If a newer buffer format is implemented, create another function, e.g.
    // DecompressV2 and keep this function for decompressing lagacy data.

    size_t bufferInOffset = 0;
    const size_t outSize = GetParameter<size_t>(bufferIn, bufferInOffset);

    m_VersionInfo = " Data is compressed using PNG Version " +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) +
                    ". Please make sure a compatible version is used for decompression.";

    png_image image;
    std::memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;

    int result = png_image_begin_read_from_memory(&image, bufferIn + bufferInOffset,
                                                  sizeIn - bufferInOffset);

    if (result == 0)
    {
        helper::Throw<std::runtime_error>("Operator", "CompressPNG", "DecompressV1",
                                          "png_image_begin_read_from_memory failed in call "
                                          "to ADIOS2 PNG Decompress." +
                                              m_VersionInfo);
    }

    // TODO might be needed from parameters?
    result = png_image_finish_read(&image, nullptr, dataOut, 0, nullptr);
    if (result == 0)
    {
        helper::Throw<std::runtime_error>("Operator", "CompressPNG", "DecompressV1",
                                          "png_image_finish_read_from_memory failed in call "
                                          "to ADIOS2 PNG Decompress." +
                                              m_VersionInfo);
    }
    return outSize;
}

void CompressPNG::CheckStatus(const int status, const std::string hint) const {}

} // end namespace compress
} // end namespace core
} // end namespace adios2
