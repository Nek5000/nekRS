/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Variable.tcc : implementation of private template functions
 *
 *  Created on: Feb 12, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_CXX11_CXX11_VARIABLE_TCC_
#define ADIOS2_BINDINGS_CXX11_CXX11_VARIABLE_TCC_

#include "Engine.h"
#include "Variable.h"
#include "adios2/core/Engine.h"
#include "adios2/helper/adiosFunctions.h"

namespace adios2
{

namespace
{

template <class T>
std::vector<typename Variable<T>::Info>
ToBlocksInfo(const std::vector<typename core::Variable<typename TypeInfo<T>::IOType>::BPInfo>
                 &coreBlocksInfo)
{
    using IOType = typename TypeInfo<T>::IOType;

    std::vector<typename Variable<T>::Info> blocksInfo;
    blocksInfo.reserve(coreBlocksInfo.size());

    for (const typename core::Variable<IOType>::BPInfo &coreBlockInfo : coreBlocksInfo)
    {
        typename Variable<T>::Info blockInfo;
        blockInfo.Start = coreBlockInfo.Start;
        blockInfo.Count = coreBlockInfo.Count;
        blockInfo.BlockID = coreBlockInfo.BlockID;
        blockInfo.Step = coreBlockInfo.Step;
        blockInfo.WriterID = coreBlockInfo.WriterID;

        blockInfo.IsReverseDims = coreBlockInfo.IsReverseDims;
        blockInfo.IsValue = coreBlockInfo.IsValue;

        if (blockInfo.IsValue)
        {
            blockInfo.Value = coreBlockInfo.Value;
        }
        else
        {
            blockInfo.Min = coreBlockInfo.Min;
            blockInfo.Max = coreBlockInfo.Max;
        }

        blocksInfo.push_back(blockInfo);
    }

    return blocksInfo;
}
} // end empty namespace

template <class T>
std::map<size_t, std::vector<typename Variable<T>::Info>>
Variable<T>::DoAllStepsBlocksInfoMap() const
{
    MinVarInfo *minBlocksInfo = nullptr;
    minBlocksInfo = m_Variable->m_Engine->MinBlocksInfo(*m_Variable, 0);
    if (!minBlocksInfo)
        throw std::logic_error("not implemented");
    std::map<size_t, std::vector<typename Variable<T>::Info>> allStepsBlocksInfo;

    size_t gotCount = 1;
    size_t curStep = 1;
    allStepsBlocksInfo.insert({0, ToBlocksInfoMin(minBlocksInfo)});
    delete (minBlocksInfo);
    while (gotCount < m_Variable->m_AvailableStepsCount)
    {
        minBlocksInfo = m_Variable->m_Engine->MinBlocksInfo(*m_Variable, curStep);
        if (minBlocksInfo)
        {
            allStepsBlocksInfo.insert({curStep, ToBlocksInfoMin(minBlocksInfo)});
            delete (minBlocksInfo);
            gotCount++;
        }
        curStep++;
    }
    return allStepsBlocksInfo;
}

template <class T>
std::vector<std::vector<typename Variable<T>::Info>> Variable<T>::DoAllStepsBlocksInfo()
{
    helper::CheckForNullptr(m_Variable, "in call to Variable<T>::AllStepsBlocksInfo");

    MinVarInfo *minBlocksInfo = nullptr;
    if (m_Variable->m_Engine)
    {
        minBlocksInfo =
            m_Variable->m_Engine->MinBlocksInfo(*m_Variable, m_Variable->m_AvailableStepsStart);
        if (minBlocksInfo)
        {
            std::vector<std::vector<typename Variable<T>::Info>> allStepsBlocksInfo;
            // PUBLIC OUTPUT
            size_t gotCount = 1;
            size_t curStep = m_Variable->m_AvailableStepsStart + 1;
            allStepsBlocksInfo.push_back(ToBlocksInfoMin(minBlocksInfo));
            delete (minBlocksInfo);
            while (gotCount < m_Variable->m_AvailableStepsCount)
            {
                minBlocksInfo = m_Variable->m_Engine->MinBlocksInfo(*m_Variable, curStep);
                if (minBlocksInfo)
                {
                    allStepsBlocksInfo.push_back(ToBlocksInfoMin(minBlocksInfo));
                    delete (minBlocksInfo);
                    gotCount++;
                }
                curStep++;
            }
            return allStepsBlocksInfo;
        }
    }

    // PRIVATE INPUT
    const std::vector<std::vector<typename core::Variable<IOType>::BPInfo>> coreAllStepsBlocksInfo =
        m_Variable->AllStepsBlocksInfo();

    // PUBLIC OUTPUT
    std::vector<std::vector<typename Variable<T>::Info>> allStepsBlocksInfo(
        coreAllStepsBlocksInfo.size());

    size_t relativeStep = 0;
    for (const auto &coreBlocksInfo : coreAllStepsBlocksInfo)
    {
        allStepsBlocksInfo[relativeStep] = ToBlocksInfo<T>(coreBlocksInfo);
        ++relativeStep;
    }
    return allStepsBlocksInfo;
}

template <typename T>
std::vector<typename Variable<T>::Info>
Variable<T>::ToBlocksInfoMin(const MinVarInfo *coreVarInfo) const
{
    auto coreBlocksInfo = coreVarInfo->BlocksInfo;
    size_t Step = coreVarInfo->Step;

    std::vector<typename Variable<T>::Info> blocksInfo;
    blocksInfo.reserve(coreBlocksInfo.size());

    for (auto &coreBlockInfo : coreBlocksInfo)
    {
        typename Variable<T>::Info blockInfo;

        blockInfo.Step = Step;
        if (coreVarInfo->Shape)
        {
            blockInfo.Start.reserve(coreVarInfo->Dims);
            blockInfo.Count.reserve(coreVarInfo->Dims);
            if (coreVarInfo->WasLocalValue)
            {
                /* Start and count are really values, not pointers */
                blockInfo.Start.push_back((size_t)coreBlockInfo.Start);
                blockInfo.Count.push_back((size_t)coreBlockInfo.Count);
            }
            else
            {
                for (int i = 0; i < coreVarInfo->Dims; i++)
                {
                    blockInfo.Start.push_back(coreBlockInfo.Start[i]);
                    blockInfo.Count.push_back(coreBlockInfo.Count[i]);
                }
            }
        }
        else
        {
            blockInfo.Count.reserve(coreVarInfo->Dims);
            for (int i = 0; i < coreVarInfo->Dims; i++)
            {
                blockInfo.Count.push_back(coreBlockInfo.Count[i]);
            }
        }
        blockInfo.WriterID = coreBlockInfo.WriterID;

        blockInfo.IsValue = coreVarInfo->IsValue;
        blockInfo.IsReverseDims = coreVarInfo->IsReverseDims;
        if (blockInfo.IsValue)
        {
            if (std::is_same<T, std::string>::value)
            {
                std::string *Tmp = (std::string *)&blockInfo.Value;
                Tmp->assign(*(const char **)coreBlockInfo.BufferP);
            }
            else
            {
                blockInfo.Value = *((T *)coreBlockInfo.BufferP);
            }
        }
        else
        {
            blockInfo.Min = *(T *)&coreBlockInfo.MinMax.MinUnion;
            blockInfo.Max = *(T *)&coreBlockInfo.MinMax.MaxUnion;
        }
        blockInfo.BlockID = coreBlockInfo.BlockID;
        blocksInfo.push_back(blockInfo);
    }

    return blocksInfo;
}

template <typename T>
std::string ToString(const Variable<T> &variable)
{
    return std::string("Variable<") + variable.Type() + ">(Name: \"" + variable.Name() + "\")";
}

#if defined(ADIOS2_HAVE_KOKKOS) || defined(ADIOS2_HAVE_GPU_SUPPORT)
template <class T>
void Variable<T>::SetArrayLayout(const ArrayOrdering layout)
{
    m_Variable->SetArrayLayout(layout);
}

template <class T>
adios2::ArrayOrdering Variable<T>::GetArrayLayout()
{
    return m_Variable->GetArrayLayout();
}
#endif

namespace detail
{
// Span
template <class T>
Span<T>::Span(core::Span<IOType> *coreSpan) : m_Span(coreSpan)
{
}

template <class T>
size_t Span<T>::size() const noexcept
{
    return m_Span->Size();
}

template <class T>
T *Span<T>::data() const noexcept
{
    return reinterpret_cast<T *>(m_Span->Data());
}

template <class T>
T &Span<T>::at(const size_t position)
{
    IOType &data = m_Span->At(position);
    return reinterpret_cast<T &>(data);
}

template <class T>
const T &Span<T>::at(const size_t position) const
{
    const IOType &data = m_Span->At(position);
    return reinterpret_cast<const T &>(data);
}

template <class T>
T &Span<T>::operator[](const size_t position)
{
    IOType &data = (*m_Span)[position];
    return reinterpret_cast<T &>(data);
}

template <class T>
const T &Span<T>::operator[](const size_t position) const
{
    const IOType &data = (*m_Span)[position];
    return reinterpret_cast<const T &>(data);
}

} // end namespace detail

} // end namespace adios2

#endif /* ADIOS2_BINDINGS_CXX11_CXX11_VARIABLE_TCC_ */
