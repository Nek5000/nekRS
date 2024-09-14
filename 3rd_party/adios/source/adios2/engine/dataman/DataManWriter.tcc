/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataManWriter.tcc
 *
 *  Created on: Jan 10, 2017
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_DATAMAN_DATAMANWRITER_TCC_
#define ADIOS2_ENGINE_DATAMAN_DATAMANWRITER_TCC_

#include "DataManWriter.h"
#include "adios2/helper/adiosSystem.h"

namespace adios2
{
namespace core
{
namespace engine
{

template <class T>
void DataManWriter::PutSyncCommon(Variable<T> &variable, const T *values)
{
    PutDeferredCommon(variable, values);
    PerformPuts();
}

template <class T>
void DataManWriter::PutDeferredCommon(Variable<T> &variable, const T *values)
{
    auto varMemSpace = variable.GetMemorySpace(values);
    variable.SetData(values);
    if (m_IO.m_ArrayOrder == ArrayOrdering::RowMajor)
    {
        m_Serializer.PutData(variable, m_Name, CurrentStep(), m_MpiRank, varMemSpace, "");
    }
    else
    {
        Dims start = variable.m_Start;
        Dims count = variable.m_Count;
        Dims shape = variable.m_Shape;
        Dims memstart = variable.m_MemoryStart;
        Dims memcount = variable.m_MemoryCount;
        std::reverse(start.begin(), start.end());
        std::reverse(count.begin(), count.end());
        std::reverse(shape.begin(), shape.end());
        std::reverse(memstart.begin(), memstart.end());
        std::reverse(memcount.begin(), memcount.end());
        m_Serializer.PutData(variable.m_Data, variable.m_Name, shape, start, count, memstart,
                             memcount, varMemSpace, m_Name, CurrentStep(), m_MpiRank, "",
                             variable.m_Operations);
    }

    if (m_MonitorActive)
    {
        m_Monitor.AddBytes(std::accumulate(variable.m_Count.begin(), variable.m_Count.end(),
                                           sizeof(T), std::multiplies<size_t>()));
    }
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_DATAMAN_DATAMANWRITER_TCC_ */
