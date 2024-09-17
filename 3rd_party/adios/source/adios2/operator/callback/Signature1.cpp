/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Signature1.cpp
 *
 *  Created on: Oct 19, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "Signature1.h"
#include "adios2/helper/adiosLog.h"

namespace adios2
{
namespace core
{
namespace callback
{

#define declare_type(T, L)                                                                         \
    Signature1::Signature1(                                                                        \
        const std::function<void(const T *, const std::string &, const std::string &,              \
                                 const std::string &, const size_t, const Dims &, const Dims &,    \
                                 const Dims &)> &function,                                         \
        const Params &parameters)                                                                  \
    : Operator("Signature1", Operator::CALLBACK_SIGNATURE1, "callback", parameters),               \
      m_Function##L(function)                                                                      \
    {                                                                                              \
    }
ADIOS2_FOREACH_STDTYPE_2ARGS(declare_type)
#undef declare_type

#define declare_type(T, L)                                                                         \
    void Signature1::RunCallback1(const T *arg1, const std::string &arg2, const std::string &arg3, \
                                  const std::string &arg4, const size_t arg5, const Dims &arg6,    \
                                  const Dims &arg7, const Dims &arg8) const                        \
    {                                                                                              \
        if (m_Function##L)                                                                         \
        {                                                                                          \
            m_Function##L(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);                         \
        }                                                                                          \
        else                                                                                       \
        {                                                                                          \
            helper::Throw<std::runtime_error>("Operator", "Signature1", "RunCallback1",            \
                                              "Signature1 with type " + std::string(#L) +          \
                                                  " callback function failed");                    \
        }                                                                                          \
    }
ADIOS2_FOREACH_STDTYPE_2ARGS(declare_type)
#undef declare_type

size_t Signature1::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                           const DataType type, char *bufferOut)
{
    return 0;
}

size_t Signature1::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    return 0;
}

bool Signature1::IsDataTypeValid(const DataType type) const { return true; }

} // end namespace callback
} // end namespace core
} // end namespace adios2
