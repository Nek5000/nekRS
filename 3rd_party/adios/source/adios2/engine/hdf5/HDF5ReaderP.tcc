/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * HDF5ReaderP.tcc
 *
 *  Created on: Oct 30, 2017
 *      Author: jgu@lbl.gov
 */

#ifndef ADIOS2_ENGINE_HDF5_HDF5FILEREADER_TCC_
#define ADIOS2_ENGINE_HDF5_HDF5FILEREADER_TCC_

#include "HDF5ReaderP.h"

namespace adios2
{
namespace core
{
namespace engine
{

template <class T>
void HDF5ReaderP::GetSyncCommon(Variable<T> &variable, T *data)
{
    // subfile info
    /*  no good way to check it is not reference to null
        if (&variable  == nullptr) {
           return;
        }
    */
    hid_t h5Type = m_H5File.GetHDF5Type<T>();
    //    UseHDFRead(variable.m_Name, data, h5Type);
    if (m_InStreamMode)
    {
        variable.m_StepsStart = m_StreamAt;
        variable.m_StepsCount = 1;
    }
    UseHDFRead(variable, data, h5Type);
}

template <class T>
void HDF5ReaderP::GetDeferredCommon(Variable<T> &variable, T *data)
{
#ifdef NEVER
    // returns immediately
    // m_HDF53Deserializer.GetDeferredVariable(variable, data);

    if (m_InStreamMode)
    {
        variable.m_StepsStart = m_StreamAt; // current step
        variable.m_StepsCount = 1;
    }
    hid_t h5Type = m_H5File.GetHDF5Type<T>();
    UseHDFRead(variable, data, h5Type);
#else
    m_DeferredStack.push_back(variable.m_Name);
    variable.SetData(data);
#endif
}

template <class T>
std::vector<typename core::Variable<T>::BPInfo>
HDF5ReaderP::BlocksInfoCommon(const core::Variable<T> &variable) const
{
    std::vector<typename core::Variable<T>::BPInfo> blocksInfo;

    typename core::Variable<T>::BPInfo blockInfo;
    blockInfo.Start = variable.m_Start;
    blockInfo.Count = variable.m_Shape;
    if (variable.m_ShapeID == ShapeID::GlobalValue || variable.m_ShapeID == ShapeID::LocalValue)
    {
        blockInfo.IsValue = true;
    }
    else
    {
        blockInfo.IsValue = false;
    }
    blocksInfo.push_back(blockInfo);

    return blocksInfo;
}

template <class T>
std::map<size_t, std::vector<typename Variable<T>::BPInfo>>
HDF5ReaderP::GetAllStepsBlocksInfo(const Variable<T> &variable) const
{
    std::map<size_t, std::vector<typename core::Variable<T>::BPInfo>> allStepsBlocksInfo;

    for (size_t step = variable.m_AvailableStepsStart; step < variable.m_AvailableStepsCount;
         ++step)
    {
        allStepsBlocksInfo[step - variable.m_AvailableStepsStart] = BlocksInfoCommon(variable);
    }
    return allStepsBlocksInfo;
}

template <class T>
std::vector<typename Variable<T>::BPInfo> HDF5ReaderP::GetBlocksInfo(const Variable<T> &variable,
                                                                     const size_t step) const
{
    return BlocksInfoCommon(variable);
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_HDF5_HDF5FILEREADER_TCC_ */
