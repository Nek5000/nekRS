/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressSirius.cpp
 *
 *  Created on: Jul 28, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "CompressSirius.h"
#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace core
{
namespace compress
{

std::unordered_map<std::string, int> CompressSirius::m_CurrentTierMap;
std::vector<std::unordered_map<std::string, std::vector<char>>> CompressSirius::m_TierBuffersMap;
int CompressSirius::m_CurrentTier;
std::vector<std::vector<char>> CompressSirius::m_TierBuffers;
int CompressSirius::m_Tiers = 0;
bool CompressSirius::m_CurrentReadFinished = false;

CompressSirius::CompressSirius(const Params &parameters)
: Operator("sirius", COMPRESS_SIRIUS, "compress", parameters)
{
    helper::GetParameter(parameters, "Tiers", m_Tiers);
    m_TierBuffersMap.resize(m_Tiers);
    m_TierBuffers.resize(m_Tiers);
}

size_t CompressSirius::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                               const DataType varType, char *bufferOut)
{
    const uint8_t bufferVersion = 1;
    size_t bufferOutOffset = 0;

    MakeCommonHeader(bufferOut, bufferOutOffset, bufferVersion);

    const size_t ndims = blockCount.size();

    // sirius V1 metadata
    PutParameter(bufferOut, bufferOutOffset, ndims);
    for (const auto &d : blockStart)
    {
        PutParameter(bufferOut, bufferOutOffset, d);
    }
    for (const auto &d : blockCount)
    {
        PutParameter(bufferOut, bufferOutOffset, d);
    }
    PutParameter(bufferOut, bufferOutOffset, varType);
    // sirius V1 metadata end

    size_t totalInputBytes = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(varType));

    // if called from Tier 0 sub-engine, then compute tier buffers and put into
    // m_TierBuffers
    size_t bytesPerTier = totalInputBytes / m_Tiers;
    if (m_CurrentTier == 0)
    {
        for (size_t i = 0; i < m_TierBuffers.size(); i++)
        {
            m_TierBuffers[i].resize(bytesPerTier);
            std::memcpy(m_TierBuffers[i].data(), dataIn + i * bytesPerTier, bytesPerTier);
        }
    }

    // for all tiers' sub-engines, copy data from m_TierBuffers to output buffer
    std::memcpy(bufferOut + bufferOutOffset, m_TierBuffers[m_CurrentTier].data(),
                m_TierBuffers[m_CurrentTier].size());

    bufferOutOffset += bytesPerTier;

    m_CurrentTier++;
    m_CurrentTier %= m_Tiers;

    return bufferOutOffset;
}

size_t CompressSirius::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
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
        // TODO: if a Version 2 sirius buffer is being implemented, put it here
        // and keep the DecompressV1 routine for backward compatibility
    }
    else
    {
        helper::Throw<std::runtime_error>("Operator", "CompressSirius", "InverseOperate",
                                          "invalid sirius buffer version");
    }

    return 0;
}

bool CompressSirius::IsDataTypeValid(const DataType type) const
{
    if (type == DataType::Float)
    {
        return true;
    }
    return false;
}

size_t CompressSirius::DecompressV1(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    // Do NOT remove even if the buffer version is updated. Data might be still
    // in lagacy formats. This function must be kept for backward compatibility.
    // If a newer buffer format is implemented, create another function, e.g.
    // DecompressV2 and keep this function for decompressing lagacy data.

    size_t bufferInOffset = 0;
    const size_t ndims = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    Dims blockStart(ndims);
    Dims blockCount(ndims);
    for (size_t i = 0; i < ndims; ++i)
    {
        blockStart[i] = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    }
    for (size_t i = 0; i < ndims; ++i)
    {
        blockCount[i] = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    }
    const DataType type = GetParameter<DataType>(bufferIn, bufferInOffset);

    const size_t outputBytes = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(type));

    std::string blockId = helper::DimsToString(blockStart) + helper::DimsToString(blockCount);

    // decompress data and copy back to m_TierBuffers
    size_t bytesPerTier = outputBytes / m_Tiers;
    auto &currentBuffer = m_TierBuffersMap[m_CurrentTierMap[blockId]][blockId];
    auto &currentTier = m_CurrentTierMap[blockId];
    currentBuffer.resize(bytesPerTier);
    std::memcpy(currentBuffer.data(), bufferIn + bufferInOffset, bytesPerTier);

    // if called from the final tier, then merge all tier buffers and copy back
    // to dataOut
    size_t accumulatedBytes = 0;

    // TODO: it currently only copies output data back when the final tier is
    // read. However, the real Sirius algorithm should instead decide when to
    // copy back decompressed data based on required acuracy level. Once it's
    // done, it should set m_CurrentReadFinished to true to inform the MHS
    // engine that the current read is finished so that it won't read the next
    // tier.
    if (currentTier == m_Tiers - 1)
    {
        for (auto &bmap : m_TierBuffersMap)
        {
            auto &b = bmap[blockId];
            std::memcpy(dataOut + accumulatedBytes, b.data(), b.size());
            accumulatedBytes += b.size();
        }
        // set m_CurrentReadFinished to true if after the current call, the
        // required acuracy is already satisfied, so that the MHS engine knows
        // it shouldn't continue reading the next tier.
        m_CurrentReadFinished = true;
    }

    // set m_CurrentReadFinished to false if the current tier does not satisfy
    // the required acuracy, so the MHS engine will read the next tier.
    m_CurrentReadFinished = false;

    currentTier++;
    if (currentTier % m_Tiers == 0)
    {
        currentTier = 0;
    }

    if (currentTier == 0)
    {
        return outputBytes;
    }
    else
    {
        return 0;
    }
}

} // end namespace compress
} // end namespace core
} // end namespace adios2
