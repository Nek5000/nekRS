/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressSirius.h
 *
 *  Created on: Jul 28, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#ifndef ADIOS2_OPERATOR_COMPRESS_COMPRESSSIRIUS_H_
#define ADIOS2_OPERATOR_COMPRESS_COMPRESSSIRIUS_H_

#include "adios2/core/Operator.h"
#include <unordered_map>

namespace adios2
{
namespace core
{
namespace compress
{

class CompressSirius : public Operator
{

public:
    CompressSirius(const Params &parameters);

    ~CompressSirius() = default;

    size_t Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                   const DataType type, char *bufferOut) final;

    size_t InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut) final;

    bool IsDataTypeValid(const DataType type) const final;

    static bool m_CurrentReadFinished;

private:
    static int m_Tiers;

    // for compress
    static std::vector<std::vector<char>> m_TierBuffers;
    static int m_CurrentTier;

    // for decompress
    static std::vector<std::unordered_map<std::string, std::vector<char>>> m_TierBuffersMap;
    static std::unordered_map<std::string, int> m_CurrentTierMap;

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
};

} // end namespace compress
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_TRANSFORM_COMPRESSION_COMPRESSSZ_H_ */
