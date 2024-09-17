/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SscReaderNaive.tcc
 *
 *  Created on: Mar 7, 2022
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_SSCREADERNAIVE_TCC_
#define ADIOS2_ENGINE_SSCREADERNAIVE_TCC_

#include "SscReaderNaive.h"
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
SscReaderNaive::BlocksInfoCommon(const Variable<T> &variable, const size_t step) const
{
    std::vector<typename Variable<T>::BPInfo> ret;
    size_t blockID = 0;
    auto it = m_BlockMap.find(variable.m_Name);
    if (it != m_BlockMap.end())
    {
        for (const auto &v : it->second)
        {
            ret.emplace_back();
            auto &b = ret.back();
            b.Start = v.start;
            b.Count = v.count;
            b.Shape = v.shape;
            b.Step = m_CurrentStep;
            b.StepsStart = m_CurrentStep;
            b.StepsCount = 1;
            b.BlockID = blockID;
            if (m_IO.m_ArrayOrder != ArrayOrdering::RowMajor)
            {
                b.IsReverseDims = true;
            }
            if (v.shapeId == ShapeID::GlobalValue || v.shapeId == ShapeID::LocalValue)
            {
                b.IsValue = true;
                std::memcpy(reinterpret_cast<char *>(&b.Value), v.value.data(), v.value.size());
            }
            ++blockID;
        }
    }
    return ret;
}

}
}
}
}

#endif // ADIOS2_ENGINE_SSCREADERNAIVE_TCC_
