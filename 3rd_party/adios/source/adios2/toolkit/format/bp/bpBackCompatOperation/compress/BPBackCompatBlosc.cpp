/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPBlosc.cpp
 *
 *  Created on: Jun 21, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "BPBackCompatBlosc.h"

#include "adios2/helper/adiosFunctions.h"

#ifdef ADIOS2_HAVE_BLOSC2
#include "adios2/operator/compress/CompressBlosc.h"
extern "C" {
#include <blosc2.h>
}
#endif

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>

namespace adios2
{
namespace format
{

void BPBackCompatBlosc::GetMetadata(const std::vector<char> &buffer, Params &info) const noexcept
{
    size_t position = 0;
    info["InputSize"] = std::to_string(helper::ReadValue<uint64_t>(buffer, position));
    info["OutputSize"] = std::to_string(helper::ReadValue<uint64_t>(buffer, position));
}

void BPBackCompatBlosc::GetData(const char *input,
                                const helper::BlockOperationInfo &blockOperationInfo,
                                char *dataOutput) const
{
#ifdef ADIOS2_HAVE_BLOSC2
    core::compress::CompressBlosc op((Params()));
    const size_t sizeOut =
        (sizeof(size_t) == 8)
            ? static_cast<size_t>(helper::StringTo<uint64_t>(
                  blockOperationInfo.Info.at("InputSize"), "when reading Blosc input size"))
            : static_cast<size_t>(helper::StringTo<uint32_t>(
                  blockOperationInfo.Info.at("InputSize"), "when reading Blosc input size"));

    Params &info = const_cast<Params &>(blockOperationInfo.Info);
    Decompress(input, blockOperationInfo.PayloadSize, dataOutput, sizeOut, info);

#else
    throw std::runtime_error("ERROR: current ADIOS2 library didn't compile "
                             "with Blosc, can't read Blosc compressed data, in call "
                             "to Get\n");
#endif
}

#ifdef ADIOS2_HAVE_BLOSC2

size_t BPBackCompatBlosc::Decompress(const void *bufferIn, const size_t sizeIn, void *dataOut,
                                     const size_t sizeOut, Params &info) const
{
    assert(sizeIn >= sizeof(DataHeader));
    const bool isChunked = reinterpret_cast<const DataHeader *>(bufferIn)->IsChunked();

    size_t decompressedSize = 0u;
    if (isChunked)
        decompressedSize = DecompressChunkedFormat(bufferIn, sizeIn, dataOut, sizeOut, info);
    else
        decompressedSize = DecompressOldFormat(bufferIn, sizeIn, dataOut, sizeOut, info);

    return decompressedSize;
}

bool BPBackCompatBlosc::IsDataTypeValid(const DataType type) const { return true; }

size_t BPBackCompatBlosc::DecompressChunkedFormat(const void *bufferIn, const size_t sizeIn,
                                                  void *dataOut, const size_t sizeOut,
                                                  Params &info) const
{
    const DataHeader *dataPtr = reinterpret_cast<const DataHeader *>(bufferIn);
    uint32_t num_chunks = dataPtr->GetNumChunks();
    size_t inputDataSize = sizeIn - sizeof(DataHeader);

    bool isCompressed = true;
    if (num_chunks == 0)
        isCompressed = false;

    size_t inputOffset = 0u;
    size_t currentOutputSize = 0u;

    const uint8_t *inputDataBuff = reinterpret_cast<const uint8_t *>(bufferIn) + sizeof(DataHeader);

    size_t uncompressedSize = sizeOut;

    if (isCompressed)
    {
        blosc2_init();
        uint8_t *outputBuff = reinterpret_cast<uint8_t *>(dataOut);

        while (inputOffset < inputDataSize)
        {
            /* move over the size of the compressed data */
            const uint8_t *in_ptr = inputDataBuff + inputOffset;

            /** read the size of the compress block from the blosc meta data
             *
             * blosc meta data format (all little endian):
             *   - 1 byte blosc format version
             *   - 1 byte blosclz format version
             *   - 1 byte flags
             *   - 1 byte typesize
             *   - 4 byte uncompressed data size
             *   - 4 byte block size
             *   - 4 byte compressed data size
             *
             * we need only the compressed size ( source address + 12 byte)
             */
            bloscSize_t max_inputDataSize = *reinterpret_cast<const bloscSize_t *>(in_ptr + 12u);

            uint8_t *out_ptr = outputBuff + currentOutputSize;

            size_t outputChunkSize = std::min<size_t>(uncompressedSize - currentOutputSize,
                                                      static_cast<size_t>(BLOSC2_MAX_BUFFERSIZE));
            bloscSize_t max_output_size = static_cast<bloscSize_t>(outputChunkSize);

            bloscSize_t decompressdSize = blosc1_decompress(in_ptr, out_ptr, max_output_size);

            if (decompressdSize > 0)
                currentOutputSize += static_cast<size_t>(decompressdSize);
            else
            {
                throw std::runtime_error(
                    "ERROR: ADIOS2 Blosc Decompress failed. Decompressed chunk "
                    "results in zero decompressed bytes.\n");
            }
            inputOffset += static_cast<size_t>(max_inputDataSize);
        }
        blosc2_destroy();
    }
    else
    {
        std::memcpy(dataOut, inputDataBuff, inputDataSize);
        currentOutputSize = inputDataSize;
        inputOffset += inputDataSize;
    }

    assert(currentOutputSize == uncompressedSize);
    assert(inputOffset == inputDataSize);

    return currentOutputSize;
}

size_t BPBackCompatBlosc::DecompressOldFormat(const void *bufferIn, const size_t sizeIn,
                                              void *dataOut, const size_t sizeOut,
                                              Params &info) const
{
    blosc2_init();
    const int decompressedSize = blosc1_decompress(bufferIn, dataOut, sizeOut);
    blosc2_destroy();
    return static_cast<size_t>(decompressedSize);
}

#endif // ADIOS2_HAVE_BLOSC2

} // end namespace format
} // end namespace adios2
