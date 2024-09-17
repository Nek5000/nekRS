/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * InlineWriter.tcc implementation of template functions with known type
 *
 *  Created on: Nov 16, 2018
 *      Author: Aron Helser aron.helser@kitware.com
 */
#ifndef ADIOS2_ENGINE_INLINEWRITER_TCC_
#define ADIOS2_ENGINE_INLINEWRITER_TCC_

#include "InlineWriter.h"

#include <iostream>

namespace adios2
{
namespace core
{
namespace engine
{

template <class T>
void InlineWriter::PutSyncCommon(Variable<T> &variable, const T *data)
{
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Writer " << m_WriterRank << "     PutSync(" << variable.m_Name
                  << ")\n";
    }

    // PutSync really shouldn't be supported for any variable, but single value
    // variables are treated as sync, even if deferred is used.
    if (variable.m_SingleValue)
    {
        PutDeferredCommon(variable, data);
    }
    else
    {
        helper::Throw<std::invalid_argument>("Engine", "InlineWriter", "PutSyncCommon",
                                             "Put Sync is not supported.");
    }
}

template <class T>
void InlineWriter::PutDeferredCommon(Variable<T> &variable, const T *data)
{
    if (m_Verbosity == 5)
    {
        std::cout << "Inline Writer " << m_WriterRank << "     PutDeferred(" << variable.m_Name
                  << ")\n";
    }

    if (m_ResetVariables)
    {
        ResetVariables();
    }
    auto &blockInfo = variable.SetBlockInfo(data, CurrentStep());
    if (variable.m_ShapeID == ShapeID::GlobalValue || variable.m_ShapeID == ShapeID::LocalValue)
    {
        blockInfo.IsValue = true;
        blockInfo.Value = blockInfo.Data[0];
    }
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_INLINEWRITER_TCC_ */
