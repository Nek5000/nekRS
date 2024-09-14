/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * NullWriter.tcc
 *
 *  Created on: 16 Apr 19
 *      Author: Chuck Atkins chuck.atkins@kitware.com
 */
#ifndef ADIOS2_ENGINE_NULLWRITER_TCC_
#define ADIOS2_ENGINE_NULLWRITER_TCC_

#include "NullWriter.h"

namespace adios2
{
namespace core
{
namespace engine
{

#define instantiate_type(T)                                                                        \
    void NullWriter::DoPut(Variable<T> &variable, typename Variable<T>::Span &span,                \
                           const bool initialize, const T &value)                                  \
    {                                                                                              \
    }

ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(instantiate_type)
#undef instantiate_type

#define instantiate_type(T)                                                                        \
    void NullWriter::DoPutSync(Variable<T> &, const T *) {}                                        \
    void NullWriter::DoPutDeferred(Variable<T> &, const T *) {}

ADIOS2_FOREACH_STDTYPE_1ARG(instantiate_type)
#undef instantiate_type

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_NULL2_NULLWRITER_TCC_ */
