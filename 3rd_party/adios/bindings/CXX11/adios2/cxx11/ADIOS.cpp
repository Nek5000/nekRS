/* Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ADIOS.cpp : public ADIOS class using PIMPL for C++11 bindings
 * Created on: Jun 4, 2018
 *     Author: William F Godoy
 */

#include "ADIOS.h"

#include "adios2/core/ADIOS.h"
#include "adios2/core/IO.h"
#include "adios2/core/VariableStruct.h"
#include "adios2/helper/adiosFunctions.h" //CheckForNullptr

namespace adios2
{
ADIOS::ADIOS(const std::string &configFile)
: m_ADIOS(std::make_shared<core::ADIOS>(configFile, "C++"))
{
}

ADIOS::ADIOS() : ADIOS("", "C++") {}

ADIOS::ADIOS(const std::string &configFile, const std::string &hostLanguage)
: m_ADIOS(std::make_shared<core::ADIOS>(configFile, hostLanguage))
{
}

ADIOS::operator bool() const noexcept { return m_ADIOS ? true : false; }

IO ADIOS::DeclareIO(const std::string name, const ArrayOrdering ArrayOrder)
{
    CheckPointer("for io name " + name + ", in call to ADIOS::DeclareIO");
    return IO(&m_ADIOS->DeclareIO(name, ArrayOrder));
}

IO ADIOS::AtIO(const std::string name)
{
    CheckPointer("for io name " + name + ", in call to ADIOS::AtIO");
    return IO(&m_ADIOS->AtIO(name));
}

void ADIOS::FlushAll()
{
    CheckPointer("in call to ADIOS::FlushAll");
    m_ADIOS->FlushAll();
}

void ADIOS::EnterComputationBlock() noexcept
{
    CheckPointer("in call to ADIOS::EnterComputationBlock()");
    m_ADIOS->EnterComputationBlock();
}

void ADIOS::ExitComputationBlock() noexcept
{
    CheckPointer("in call to ADIOS::ExitComputationBlock()");
    m_ADIOS->ExitComputationBlock();
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
    if (op == nullptr)
    {
        return Operator("", nullptr);
    }
    else
    {
        return Operator(op->first, &op->second);
    }
}

bool ADIOS::RemoveIO(const std::string name)
{
    CheckPointer("for io name " + name + ", in call to ADIOS::RemoveIO");
    return m_ADIOS->RemoveIO(name);
}

void ADIOS::RemoveAllIOs() noexcept
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

} // end namespace adios2
