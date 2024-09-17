/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Variable.tcc : implement long template functions
 *
 *  Created on: Jan 31, 2019
 *      Author: William F Godoy
 */

#ifndef ADIOS2_CORE_VARIABLE_TCC_
#define ADIOS2_CORE_VARIABLE_TCC_

#include "Variable.h"

#include "adios2/core/Engine.h"
#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace core
{

template <class T>
Dims Variable<T>::DoCount() const
{
    auto lf_Step = [&]() -> size_t {
        auto itStep = std::next(m_AvailableStepBlockIndexOffsets.begin(), m_StepsStart);
        if (itStep == m_AvailableStepBlockIndexOffsets.end())
        {
            auto it = m_AvailableStepBlockIndexOffsets.rbegin();
            helper::Throw<std::invalid_argument>(
                "Core", "Variable", "DoCount",
                "current relative step start for variable " + m_Name +
                    " is outside the scope of available steps " + std::to_string(it->first - 1) +
                    " in call to Count");
        }
        return itStep->first - 1;
    };

    if (m_Engine != nullptr && m_SelectionType == SelectionType::WriteBlock)
    {
        auto MVI = m_Engine->MinBlocksInfo(*this, m_StepsStart);
        if (MVI)
        {
            if (m_BlockID >= MVI->BlocksInfo.size())
            {
                helper::Throw<std::invalid_argument>(
                    "Core", "Variable", "DoCount",
                    "blockID " + std::to_string(m_BlockID) +
                        " from SetBlockSelection is out of bounds for "
                        "available "
                        "blocks size " +
                        std::to_string(MVI->BlocksInfo.size()) + " for variable " + m_Name +
                        " for step " + std::to_string(m_StepsStart) +
                        ", in call to Variable<T>::Count()");
            }

            if (!MVI->WasLocalValue)
            {
                size_t *DimsPtr = (MVI->BlocksInfo)[m_BlockID].Count;
                Dims D;
                D.resize(MVI->Dims);
                for (int i = 0; i < MVI->Dims; i++)
                {
                    D[i] = DimsPtr[i];
                }
                delete MVI;
                return D;
            }
            else
            {
                delete MVI;
                return {1};
            }
        }

        const size_t step = !m_FirstStreamingStep ? m_Engine->CurrentStep() : lf_Step();

        const std::vector<typename Variable<T>::BPInfo> blocksInfo =
            m_Engine->BlocksInfo<T>(*this, step);

        if (m_BlockID >= blocksInfo.size())
        {
            helper::Throw<std::invalid_argument>(
                "Core", "Variable", "DoCount",
                "blockID " + std::to_string(m_BlockID) +
                    " from SetBlockSelection is out of bounds for available "
                    "blocks size " +
                    std::to_string(blocksInfo.size()) + " for variable " + m_Name + " for step " +
                    std::to_string(step) + ", in call to Variable<T>::Count()");
        }

        return blocksInfo[m_BlockID].Count;
    }
    return m_Count;
}

template <class T>
size_t Variable<T>::DoSelectionSize() const
{
    return helper::GetTotalSize(DoCount()) * m_StepsCount;
}

template <class T>
std::pair<T, T> Variable<T>::DoMinMax(const size_t step) const
{
    CheckRandomAccess(step, "MinMax");

    std::pair<T, T> minMax;
    minMax.first = {};
    minMax.second = {};

    if (m_Engine != nullptr)
    {

        MinMaxStruct MM;
        if (m_Engine->VariableMinMax(*this, step, MM))
        {
            if (std::is_same<T, std::string>::value)
            {
                return minMax;
            }
            else
            {
                minMax.first = *(T *)&MM.MinUnion;
                minMax.second = *(T *)&MM.MaxUnion;
                return minMax;
            }
        }
    }
    if (m_Engine != nullptr && !m_FirstStreamingStep)
    {
        const size_t stepInput = (step == DefaultSizeT) ? m_Engine->CurrentStep() : step;

        const std::vector<typename Variable<T>::BPInfo> blocksInfo =
            m_Engine->BlocksInfo<T>(*this, stepInput);

        if (blocksInfo.size() == 0)
        {
            return minMax;
        }

        if (m_ShapeID == ShapeID::LocalArray)
        {
            if (m_BlockID >= blocksInfo.size())
            {
                helper::Throw<std::invalid_argument>(
                    "Core", "Variable", "DoMinMax",
                    "BlockID " + std::to_string(m_BlockID) +
                        " does not exist for LocalArray variable " + m_Name +
                        ", in call to MinMax, Min or Maxn");
            }
            minMax.first = blocksInfo[m_BlockID].Min;
            minMax.second = blocksInfo[m_BlockID].Max;
            return minMax;
        }

        const bool isValue = ((blocksInfo.front().Shape.size() == 1 &&
                               blocksInfo.front().Shape.front() == LocalValueDim) ||
                              m_ShapeID == ShapeID::GlobalValue)
                                 ? true
                                 : false;

        minMax.first = isValue ? blocksInfo.front().Value : blocksInfo.front().Min;
        minMax.second = isValue ? blocksInfo.front().Value : blocksInfo.front().Max;

        for (const typename Variable<T>::BPInfo &blockInfo : blocksInfo)
        {
            const T minValue = isValue ? blockInfo.Value : blockInfo.Min;

            if (helper::LessThan<T>(minValue, minMax.first))
            {
                minMax.first = minValue;
            }

            const T maxValue = isValue ? blockInfo.Value : blockInfo.Max;

            if (helper::GreaterThan<T>(maxValue, minMax.second))
            {
                minMax.second = maxValue;
            }
        }
    }
    else
    {
        minMax.first = m_Min;
        minMax.second = m_Max;
    }
    return minMax;
}

template <class T>
std::vector<std::vector<typename Variable<T>::BPInfo>> Variable<T>::DoAllStepsBlocksInfo() const
{
    if (m_Engine == nullptr)
    {
        helper::Throw<std::invalid_argument>("Core", "Variable", "DoAllStepsBlocksInfo",
                                             "from variable " + m_Name +
                                                 " function is only valid in read mode, in "
                                                 "call to Variable<T>::AllBlocksInfo");
    }

    if (!m_FirstStreamingStep)
    {
        helper::Throw<std::invalid_argument>("Core", "Variable", "DoAllStepsBlocksInfo",
                                             "from variable " + m_Name +
                                                 " function is not valid in "
                                                 "random-access read mode "
                                                 "(BeginStep/EndStep), in "
                                                 "call to Variable<T>::AllBlocksInfo");
    }

    return m_Engine->AllRelativeStepsBlocksInfo(*this);
}

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_VARIABLE_TCC_ */
