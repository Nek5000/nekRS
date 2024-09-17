/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * NullReader.tcc
 *
 */
#ifndef ADIOS2_ENGINE_NULLREADER_TCC_
#define ADIOS2_ENGINE_NULLREADER_TCC_

#include "NullReader.h"

namespace adios2
{
namespace core
{
namespace engine
{

template <>
inline void NullReader::GetSyncCommon(Variable<std::string> &variable, std::string *data)
{
    variable.m_Data = data;
}

template <class T>
inline void NullReader::GetSyncCommon(Variable<T> &variable, T *data)
{
    variable.m_Data = data;
}

template <class T>
void NullReader::GetDeferredCommon(Variable<T> &variable, T *data)
{
    variable.m_Data = data;
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_NULLREADER_TCC_ */
