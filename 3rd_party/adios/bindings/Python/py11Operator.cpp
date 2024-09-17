/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11Operator.cpp
 *
 *  Created on: Dec 12, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "py11Operator.h"

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace py11
{

Operator::operator bool() const noexcept { return m_Parameters != nullptr; }

std::string Operator::Type() const noexcept
{
    if (m_Parameters == nullptr)
    {
        helper::Log("PythonAPI", "Operator", "Type()", "Operator is nullptr",
                    helper::LogMode::EXCEPTION);
    }
    return m_Type;
}

void Operator::SetParameter(const std::string key, const std::string value)
{
    if (m_Parameters == nullptr)
    {
        helper::Log("PythonAPI", "Operator", "SetParameter()", "Operator is nullptr",
                    helper::LogMode::EXCEPTION);
    }
    (*m_Parameters)[key] = value;
}

Params &Operator::Parameters() const
{
    if (m_Parameters == nullptr)
    {
        helper::Log("PythonAPI", "Operator", "Parameter()", "Operator is nullptr",
                    helper::LogMode::EXCEPTION);
    }
    return *m_Parameters;
}

// PRIVATE
Operator::Operator(const std::string &type, Params *params) : m_Parameters(params), m_Type(type) {}

} // end namespace py11
} // end namespace adios2
