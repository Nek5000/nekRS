/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11ADIOS.cpp
 *
 *  Created on: Mar 13, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "py11ADIOS.h"

namespace adios2
{
namespace py11
{

ADIOS::ADIOS(const std::string &configFile)
: m_ADIOS(std::make_shared<adios2::core::ADIOS>(configFile, "Python"))
{
}

ADIOS::ADIOS() : m_ADIOS(std::make_shared<adios2::core::ADIOS>("", "Python")) {}

ADIOS::operator bool() const noexcept { return m_ADIOS ? true : false; }

IO ADIOS::DeclareIO(const std::string name)
{
    CheckPointer("for io name " + name + ", in call to ADIOS::DeclareIO");
    return IO(&m_ADIOS->DeclareIO(name));
}

IO ADIOS::AtIO(const std::string name)
{
    CheckPointer("for io name " + name + ", in call to ADIOS::AtIO");
    return IO(&m_ADIOS->AtIO(name));
}

Operator ADIOS::DefineOperator(const std::string name, const std::string type,
                               const Params &parameters)
{
    CheckPointer("for operator name " + name + ", in call to ADIOS::DefineOperator");
    auto op = &m_ADIOS->DefineOperator(name, type, parameters);
    return Operator(op->first, &op->second);
}

Operator ADIOS::InquireOperator(const std::string name)
{
    CheckPointer("for operator name " + name + ", in call to InquireOperator");
    auto op = m_ADIOS->InquireOperator(name);
    return Operator(op->first, &op->second);
}

void ADIOS::FlushAll()
{
    CheckPointer("in call to ADIOS::FlushAll");
    m_ADIOS->FlushAll();
}

bool ADIOS::RemoveIO(const std::string name)
{
    CheckPointer("in call to ADIOS::RemoveIO");
    return m_ADIOS->RemoveIO(name);
}

void ADIOS::RemoveAllIOs()
{
    CheckPointer("in call to ADIOS::RemoveAllIOs");
    m_ADIOS->RemoveAllIOs();
}

// PRIVATE
void ADIOS::CheckPointer(const std::string hint)
{
    if (!m_ADIOS)
    {
        throw std::invalid_argument("ERROR: invalid ADIOS object, did you call "
                                    "any of the ADIOS explicit "
                                    "constructors?, " +
                                    hint + "\n");
    }
}

} // end namespace py11
} // end namespace adios2
