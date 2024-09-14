/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP3Reader.tcc
 *
 *  Created on: Feb 27, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_ENGINE_BP3_BP3READER_TCC_
#define ADIOS2_ENGINE_BP3_BP3READER_TCC_

#include "BP3Reader.h"

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace core
{
namespace engine
{

template <>
inline void BP3Reader::GetSyncCommon(Variable<std::string> &variable, std::string *data)
{
    m_BP3Deserializer.GetValueFromMetadata(variable, data);
}

template <class T>
inline void BP3Reader::GetSyncCommon(Variable<T> &variable, T *data)
{
    if (variable.m_SingleValue)
    {
        m_BP3Deserializer.GetValueFromMetadata(variable, data);
        return;
    }

    typename Variable<T>::BPInfo &blockInfo =
        m_BP3Deserializer.InitVariableBlockInfo(variable, data);
    m_BP3Deserializer.SetVariableBlockInfo(variable, blockInfo);
    ReadVariableBlocks(variable);
    variable.m_BlocksInfo.pop_back();
}

template <class T>
void BP3Reader::GetDeferredCommon(Variable<T> &variable, T *data)
{
    // cheap
    if (variable.m_SingleValue)
    {
        m_BP3Deserializer.GetValueFromMetadata(variable, data);
        return;
    }

    // returns immediately without populating data
    m_BP3Deserializer.InitVariableBlockInfo(variable, data);
    m_BP3Deserializer.m_DeferredVariables.insert(variable.m_Name);
}

template <class T>
void BP3Reader::ReadVariableBlocks(Variable<T> &variable)
{
    const bool profile = m_BP3Deserializer.m_Profiler.m_IsActive;

    for (typename Variable<T>::BPInfo &blockInfo : variable.m_BlocksInfo)
    {
        T *originalBlockData = blockInfo.Data;

        for (const auto &stepPair : blockInfo.StepBlockSubStreamsInfo)
        {
            for (const helper::SubStreamBoxInfo &subStreamBoxInfo : stepPair.second)
            {
                if (subStreamBoxInfo.ZeroBlock)
                {
                    continue;
                }

                // check if subfile is already opened
                if (m_SubFileManager.m_Transports.count(subStreamBoxInfo.SubStreamID) == 0)
                {
                    const std::string subFileName = m_BP3Deserializer.GetBPSubFileName(
                        m_Name, subStreamBoxInfo.SubStreamID,
                        m_BP3Deserializer.m_Minifooter.HasSubFiles, true);

                    m_SubFileManager.OpenFileID(subFileName, subStreamBoxInfo.SubStreamID,
                                                Mode::Read, {{"transport", "File"}}, profile);
                }

                char *buffer = nullptr;
                size_t payloadSize = 0, payloadStart = 0;

                m_BP3Deserializer.PreDataRead(variable, blockInfo, subStreamBoxInfo, buffer,
                                              payloadSize, payloadStart, 0);

                m_SubFileManager.ReadFile(buffer, payloadSize, payloadStart,
                                          subStreamBoxInfo.SubStreamID);

                m_BP3Deserializer.PostDataRead(variable, blockInfo, subStreamBoxInfo,
                                               m_IO.m_ArrayOrder == ArrayOrdering::RowMajor, 0);
            } // substreams loop
            // advance pointer to next step
            blockInfo.Data += helper::GetTotalSize(blockInfo.Count);
        } // steps loop
        blockInfo.Data = originalBlockData;
    } // deferred blocks loop
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_BP3_BP3READER_TCC_ */
