/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DaosReader.tcc
 *
 */

#ifndef ADIOS2_ENGINE_DAOS_DAOSREADER_TCC_
#define ADIOS2_ENGINE_DAOS_DAOSREADER_TCC_

#include "DaosReader.h"

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace core
{
namespace engine
{

inline void DaosReader::GetSyncCommon(VariableBase &variable, void *data)
{
    bool need_sync = m_BP5Deserializer->QueueGet(variable, data);
    if (need_sync)
        PerformGets();
}

void DaosReader::GetDeferredCommon(VariableBase &variable, void *data)
{
    (void)m_BP5Deserializer->QueueGet(variable, data);
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_DAOS_DAOSREADER_TCC_ */
