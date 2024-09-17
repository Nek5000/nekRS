/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscReaderGeneric.tcc
 *
 *  Created on: Mar 3, 2022
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCREADERGENERIC_TCC_
#define ADIOS2_ENGINE_SSCREADERGENERIC_TCC_

#include "SscReaderGeneric.h"
#include "adios2/helper/adiosMemory.h"

namespace adios2
{
namespace core
{
namespace engine
{
namespace ssc
{

template <typename T>
std::vector<typename Variable<T>::BPInfo>
SscReaderGeneric::BlocksInfoCommon(const Variable<T> &variable, const size_t step) const
{
    std::vector<typename Variable<T>::BPInfo> ret;
    size_t blockID = 0;
    for (int i = 0; i < static_cast<int>(m_GlobalWritePattern.size()); ++i)
    {
        for (auto &v : m_GlobalWritePattern[i])
        {
            if (v.name == variable.m_Name)
            {
                ret.emplace_back();
                auto &b = ret.back();
                b.Start = v.start;
                b.Count = v.count;
                b.Shape = v.shape;
                b.Step = m_CurrentStep;
                b.StepsStart = m_CurrentStep;
                b.StepsCount = 1;
                b.WriterID = i;
                b.BlockID = blockID;
                if (m_IO.m_ArrayOrder != ArrayOrdering::RowMajor)
                {
                    b.IsReverseDims = true;
                }
                if (v.shapeId == ShapeID::GlobalValue || v.shapeId == ShapeID::LocalValue)
                {
                    b.IsValue = true;
                    if (m_CurrentStep == 0 || m_WriterDefinitionsLocked == false ||
                        m_ReaderSelectionsLocked == false)
                    {
                        std::memcpy(reinterpret_cast<char *>(&b.Value), v.value.data(),
                                    v.value.size());
                    }
                    else
                    {
                        std::memcpy(reinterpret_cast<char *>(&b.Value),
                                    m_Buffer.data() + v.bufferStart, v.bufferCount);
                    }
                }
                ++blockID;
            }
        }
    }
    return ret;
}

}
}
}
}

#endif // ADIOS2_ENGINE_SSCREADERGENERIC_TCC_
