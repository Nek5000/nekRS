/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressBlosc.h
 *
 *  Created on: Jun 18, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 *              Rene Widera r.widera@hzdr.de
 */

#ifndef ADIOS2_OPERATOR_COMPRESS_COMPRESSBLOSC_H_
#define ADIOS2_OPERATOR_COMPRESS_COMPRESSBLOSC_H_

#include <map>
#include <set>

#include "adios2/core/Operator.h"

#if defined(_MSC_VER)
#define ADIOS2_CLASS_PACKED(name) __pragma(pack(push, 1)) class name
#define ADIOS2_CLASS_PACKED_SUFFIX __pragma(pack(pop))
#else
#define ADIOS2_CLASS_PACKED(name) class __attribute__((packed)) name
#define ADIOS2_CLASS_PACKED_SUFFIX
#endif

namespace adios2
{
namespace core
{
namespace compress
{

class CompressBlosc : public Operator
{

public:
    /**
     * Unique constructor
     */
    CompressBlosc(const Params &parameters);

    ~CompressBlosc() = default;

    /**
     * @param dataIn
     * @param blockStart
     * @param blockCount
     * @param type
     * @param bufferOut format will be: 'DataHeader ; (BloscCompressedChunk |
     * UncompressedData), [ BloscCompressedChunk, ...]'
     * @return size of compressed buffer in bytes
     */
    size_t Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                   const DataType type, char *bufferOut) final;

    /**
     * @param bufferIn
     * @param sizeIn
     * @param dataOut
     * @return size of decompressed buffer
     */
    size_t InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut) final;

    bool IsDataTypeValid(const DataType type) const final;

    size_t GetHeaderSize() const;

private:
    using bloscSize_t = int32_t;
    size_t headerSize = 0;
    size_t bufferDecompressedSize = 0;

    /** Decompress chunked data */
    size_t DecompressChunkedFormat(const char *bufferIn, const size_t sizeIn, char *dataOut,
                                   const size_t sizeOut);

    /** Decompress data written before ADIOS2 supported large variables larger
     * 2GiB. */
    size_t DecompressOldFormat(const char *bufferIn, const size_t sizeIn, char *dataOut,
                               const size_t sizeOut) const;

    /**
     * Decompress function for V1 buffer. Do NOT remove even if the buffer
     * version is updated. Data might be still in lagacy formats. This function
     * must be kept for backward compatibility
     * @param bufferIn : compressed data buffer (V1 only)
     * @param sizeIn : number of bytes in bufferIn
     * @param dataOut : decompressed data buffer
     * @return : number of bytes in dataOut
     */
    size_t DecompressV1(const char *bufferIn, const size_t sizeIn, char *dataOut);

    ADIOS2_CLASS_PACKED(DataHeader)
    {
        /** compatible to the first 4 byte of blosc2 header
         *
         *   blosc2 meta data format (all little endian):
         *   - 1 byte blosc2 format version
         *   - 1 byte blosc2lz format version
         *   - 1 byte flags
         *   - 1 byte typesize
         *
         * If zero we writing the new adios blosc2 format which can handle more
         * than 2GiB data chunks.
         */
        uint32_t format = 0u;
        /** number of blosc2 chunks within the data blob
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

    static const std::map<std::string, uint32_t> m_Shuffles;
    static const std::set<std::string> m_Compressors;

    std::string m_VersionInfo;
};

} // end namespace compress
} // end namespace core
} // end namespace adios2

#undef ADIOS2_CLASS_PACKED
#undef ADIOS2_CLASS_PACKED_SUFFIX

#endif /* ADIOS2_OPERATOR_COMPRESS_COMPRESSBLOSC_H_ */
