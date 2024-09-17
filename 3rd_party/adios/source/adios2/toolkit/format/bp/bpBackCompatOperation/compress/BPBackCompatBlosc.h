/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPBlosc.h
 *
 *  Created on: Jun 21, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPOPERATION_COMPRESS_BPBLOSC_H_
#define ADIOS2_TOOLKIT_FORMAT_BP_BPOPERATION_COMPRESS_BPBLOSC_H_

#include "adios2/toolkit/format/bp/bpBackCompatOperation/BPBackCompatOperation.h"

#if defined(_MSC_VER)
#define ADIOS2_CLASS_PACKED(name) __pragma(pack(push, 1)) class name
#define ADIOS2_CLASS_PACKED_SUFFIX __pragma(pack(pop))
#else
#define ADIOS2_CLASS_PACKED(name) class __attribute__((packed)) name
#define ADIOS2_CLASS_PACKED_SUFFIX
#endif

namespace adios2
{
namespace format
{

class BPBackCompatBlosc : public BPBackCompatOperation
{
public:
    BPBackCompatBlosc() = default;

    ~BPBackCompatBlosc() = default;

    void GetMetadata(const std::vector<char> &buffer, Params &info) const noexcept final;

    void GetData(const char *input, const helper::BlockOperationInfo &blockOperationInfo,
                 char *dataOutput) const final;

private:
    /**
     * Decompression signature for legacy libraries that use void*
     * @param bufferIn
     * @param sizeIn
     * @param dataOut
     * @param dimensions
     * @param type
     * @return size of decompressed buffer in bytes
     */
    size_t Decompress(const void *bufferIn, const size_t sizeIn, void *dataOut,
                      const size_t sizeOut, Params &info) const;

    bool IsDataTypeValid(const DataType type) const;

    using bloscSize_t = int32_t;

    /** Decompress chunked data */
    size_t DecompressChunkedFormat(const void *bufferIn, const size_t sizeIn, void *dataOut,
                                   const size_t sizeOut, Params &info) const;

    /** Decompress data written before ADIOS2 supported large variables larger
     * 2GiB. */
    size_t DecompressOldFormat(const void *bufferIn, const size_t sizeIn, void *dataOut,
                               const size_t sizeOut, Params &info) const;

    ADIOS2_CLASS_PACKED(DataHeader)
    {
        /** compatible to the first 4 byte of blosc header
         *
         *   blosc meta data format (all little endian):
         *   - 1 byte blosc format version
         *   - 1 byte blosclz format version
         *   - 1 byte flags
         *   - 1 byte typesize
         *
         * If zero we writing the new adios blosc format which can handle more
         * than 2GiB data chunks.
         */
        uint32_t format = 0u;
        /** number of blosc chunks within the data blob
         *
         * If zero the data is not compressed and must be decompressed by using
         * 'memcpy'
         */
        uint32_t numberOfChunks = 0u;

    public:
        void SetNumChunks(const uint32_t numChunks) { numberOfChunks = numChunks; }
        uint32_t GetNumChunks() const { return numberOfChunks; }

        bool IsChunked() const { return format == 0; }
    }
    ADIOS2_CLASS_PACKED_SUFFIX;
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BP_BPOPERATION_COMPRESS_BPBLOSC_H_ */
