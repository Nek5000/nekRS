/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Operator.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "Operator.h"
#include "adios2/helper/adiosFunctions.h"

#include <iostream>

namespace adios2
{
namespace core
{

Operator::Operator(const std::string &typeString, const OperatorType typeEnum,
                   const std::string &category, const Params &parameters)
: m_TypeString(typeString), m_TypeEnum(typeEnum), m_Category(category),
  m_Parameters(helper::LowerCaseParams(parameters))
{
}

void Operator::SetParameter(const std::string key, const std::string value) noexcept
{
    m_Parameters[helper::LowerCase(key)] = value;
}

Params &Operator::GetParameters() noexcept { return m_Parameters; }

void Operator::SetAccuracy(const adios2::Accuracy &a) noexcept { m_AccuracyRequested = a; }
adios2::Accuracy Operator::GetAccuracy() const noexcept { return m_AccuracyProvided; }

#define declare_type(T)                                                                            \
                                                                                                   \
    void Operator::RunCallback1(const T *arg0, const std::string &arg1, const std::string &arg2,   \
                                const std::string &arg3, const size_t arg4, const Dims &arg5,      \
                                const Dims &arg6, const Dims &arg7) const                          \
    {                                                                                              \
        CheckCallbackType("Callback1");                                                            \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void Operator::RunCallback2(void *arg0, const std::string &arg1, const std::string &arg2,
                            const std::string &arg3, const size_t arg4, const Dims &arg5,
                            const Dims &arg6, const Dims &arg7) const
{
    CheckCallbackType("Callback2");
}

// PROTECTED

Dims Operator::ConvertDims(const Dims &dimensions, const DataType type, const size_t targetDims,
                           const bool enforceDims, const size_t defaultDimSize) const
{

    if (targetDims < 1)
    {
        helper::Throw<std::invalid_argument>("Core", "Operator", "ConvertDims",
                                             "only accepts targetDims > 0");
    }

    Dims ret = dimensions;

    while (true)
    {
        auto it = std::find(ret.begin(), ret.end(), 1);
        if (it == ret.end())
        {
            break;
        }
        else
        {
            ret.erase(it);
        }
    }

    while (ret.size() > targetDims)
    {
        ret[1] *= ret[0];
        ret.erase(ret.begin());
    }

    while (enforceDims && ret.size() < targetDims)
    {
        ret.insert(ret.begin(), defaultDimSize);
    }

    if (type == helper::GetDataType<std::complex<float>>() ||
        type == helper::GetDataType<std::complex<double>>())
    {
        ret.back() *= 2;
    }
    return ret;
}

size_t Operator::GetHeaderSize() const { return 0; }

size_t Operator::GetEstimatedSize(const size_t ElemCount, const size_t ElemSize, const size_t ndims,
                                  const size_t *dims) const
{
    return ElemCount * ElemSize + 128;
};

// PRIVATE
void Operator::CheckCallbackType(const std::string type) const
{
    if (m_TypeString != type)
    {
        helper::Throw<std::invalid_argument>("Core", "Operator", "CheckCallbackType",
                                             "operator of type " + m_TypeString +
                                                 " doesn't match expected callback type " + type +
                                                 " arguments");
    }
}

} // end namespace core
} // end namespace adios2
