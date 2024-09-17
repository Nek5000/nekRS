/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SkeletonWriter.tcc implementation of template functions with known type
 *
 *  Created on: Jan 04, 2018
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */
#ifndef ADIOS2_ENGINE_SKELETONWRITER_TCC_
#define ADIOS2_ENGINE_SKELETONWRITER_TCC_

#include "SkeletonWriter.h"

#include <iostream>

namespace adios2
{
namespace core
{
namespace engine
{

template <class T>
void SkeletonWriter::PutSyncCommon(Variable<T> &variable,
                                   const typename Variable<T>::BPInfo &blockInfo)
{
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Writer " << m_WriterRank << "     PutSync(" << variable.m_Name
                  << ")\n";
    }
}

template <class T>
void SkeletonWriter::PutDeferredCommon(Variable<T> &variable, const T *data)
{
    variable.SetBlockInfo(data, CurrentStep());

    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Writer " << m_WriterRank << "     PutDeferred(" << variable.m_Name
                  << ")\n";
    }
    m_NeedPerformPuts = true;
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_SKELETONWRITER_TCC_ */
