/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CampaignReader.tcc
 *
 *  Created on: May 15, 2023
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#ifndef ADIOS2_ENGINE_CAMPAIGNREADER_TCC_
#define ADIOS2_ENGINE_CAMPAIGNREADER_TCC_

#include "CampaignReader.h"

#include <iostream>

namespace adios2
{
namespace core
{
namespace engine
{

template <class T>
inline Variable<T> CampaignReader::DuplicateVariable(Variable<T> *variable, IO &io,
                                                     std::string &name, VarInternalInfo &vii)
{
    auto &v = io.DefineVariable<T>(name, variable->Shape());
    v.m_AvailableStepsCount = variable->GetAvailableStepsCount();
    v.m_AvailableStepsStart = variable->GetAvailableStepsStart();
    v.m_ShapeID = variable->m_ShapeID;
    v.m_SingleValue = variable->m_SingleValue;
    v.m_ReadAsJoined = variable->m_ReadAsJoined;
    v.m_ReadAsLocalValue = variable->m_ReadAsLocalValue;
    v.m_RandomAccess = variable->m_RandomAccess;
    v.m_MemSpace = variable->m_MemSpace;
    v.m_JoinedDimPos = variable->m_JoinedDimPos;
    v.m_AvailableStepBlockIndexOffsets = variable->m_AvailableStepBlockIndexOffsets;
    v.m_AvailableShapes = variable->m_AvailableShapes;
    v.m_Min = variable->m_Min;
    v.m_Max = variable->m_Max;
    v.m_Value = variable->m_Value;
    v.m_StepsStart = variable->m_StepsStart;
    v.m_StepsCount = variable->m_StepsCount;
    v.m_Start = variable->m_Start;
    v.m_Count = variable->m_Count;

    v.m_Engine = this; // Variable::Shape() uses this member to call engine
    vii.originalVar = static_cast<void *>(variable);
    m_VarInternalInfo.emplace(name, vii);
    return v;
}

template <class T>
inline Attribute<T> CampaignReader::DuplicateAttribute(Attribute<T> *attribute, IO &io,
                                                       std::string &name)
{
    if (attribute->m_IsSingleValue)
    {
        auto &a = io.DefineAttribute<T>(name, attribute->m_DataSingleValue);
        return a;
    }
    auto &a = io.DefineAttribute<T>(name, attribute->m_DataArray.data(), attribute->m_Elements);
    return a;
}

template <class T>
inline std::pair<Variable<T> *, Engine *>
CampaignReader::TranslateToActualVariable(Variable<T> &variable)
{
    auto it = m_VarInternalInfo.find(variable.m_Name);
    Variable<T> *v = reinterpret_cast<Variable<T> *>(it->second.originalVar);
    Engine *e = m_Engines[it->second.engineIdx];
    v->m_SelectionType = variable.m_SelectionType;
    v->m_Start = variable.m_Start;
    v->m_Count = variable.m_Count;
    v->m_StepsStart = variable.m_StepsStart;
    v->m_StepsCount = variable.m_StepsCount;
    v->m_BlockID = variable.m_BlockID;
    v->m_MemoryStart = variable.m_MemoryStart;
    v->m_MemoryCount = variable.m_MemoryCount;
    v->m_MemSpace = variable.m_MemSpace;
    return std::make_pair(v, e);
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif // ADIOS2_ENGINE_CAMPAIGNREADER_TCC_
