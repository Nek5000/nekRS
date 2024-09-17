/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO.cpp :
 *
 *  Created on: Jun 4, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "IO.h"
#include "IO.tcc"

#include "adios2/core/IO.h"

namespace adios2
{

IO::operator bool() const noexcept { return (m_IO == nullptr) ? false : true; }

std::string IO::Name() const
{
    helper::CheckForNullptr(m_IO, "in call to IO::InConfigFile");
    return m_IO->m_Name;
}

bool IO::InConfigFile() const
{
    helper::CheckForNullptr(m_IO, "in call to IO::InConfigFile");
    return m_IO->InConfigFile();
}

void IO::SetEngine(const std::string engineType)
{
    helper::CheckForNullptr(m_IO, "in call to IO::SetEngine");
    m_IO->SetEngine(engineType);
}

void IO::SetParameter(const std::string key, const std::string value)
{
    helper::CheckForNullptr(m_IO, "in call to IO::SetParameter");
    m_IO->SetParameter(key, value);
}

void IO::SetParameters(const Params &parameters)
{
    helper::CheckForNullptr(m_IO, "in call to IO::SetParameters");
    m_IO->SetParameters(parameters);
}

void IO::SetParameters(const std::string &parameters)
{
    helper::CheckForNullptr(m_IO, "in call to IO::SetParameters");
    m_IO->SetParameters(parameters);
}

void IO::ClearParameters()
{
    helper::CheckForNullptr(m_IO, "in call to IO::ClearParameters");
    m_IO->ClearParameters();
}

Params IO::Parameters() const
{
    helper::CheckForNullptr(m_IO, "in call to IO:::Parameters");
    return m_IO->m_Parameters;
}

size_t IO::AddTransport(const std::string type, const Params &parameters)
{
    helper::CheckForNullptr(m_IO, "in call to IO::AddTransport");
    return m_IO->AddTransport(type, parameters);
}

void IO::SetTransportParameter(const size_t transportIndex, const std::string key,
                               const std::string value)
{
    helper::CheckForNullptr(m_IO, "in call to IO::SetTransportParameter");
    m_IO->SetTransportParameter(transportIndex, key, value);
}

bool IO::RemoveVariable(const std::string &name)
{
    helper::CheckForNullptr(m_IO, "in call to IO::RemoveVariable");
    return m_IO->RemoveVariable(name);
}

void IO::RemoveAllVariables()
{
    helper::CheckForNullptr(m_IO, "in call to IO::RemoveAllVariables");
    m_IO->RemoveAllVariables();
}

bool IO::RemoveAttribute(const std::string &name)
{
    helper::CheckForNullptr(m_IO, "in call to IO::RemoveAttribute");
    return m_IO->RemoveAttribute(name);
}

void IO::RemoveAllAttributes()
{
    helper::CheckForNullptr(m_IO, "in call to IO::RemoveAllAttributes");
    m_IO->RemoveAllAttributes();
}

Engine IO::Open(const std::string &name, const Mode mode)
{
    helper::CheckForNullptr(m_IO, "for engine " + name + ", in call to IO::Open");
    return Engine(&m_IO->Open(name, mode));
}
Group IO::InquireGroup(char delimiter) { return Group(&m_IO->CreateGroup(delimiter)); };
void IO::FlushAll()
{
    helper::CheckForNullptr(m_IO, "in call to IO::FlushAll");
    m_IO->FlushAll();
}

std::map<std::string, Params> IO::AvailableVariables(bool namesOnly)
{
    helper::CheckForNullptr(m_IO, "in call to IO::AvailableVariables");
    if (namesOnly)
    {
        return m_IO->GetAvailableVariables({"name"});
    }
    else
    {
        return m_IO->GetAvailableVariables();
    }
}

std::map<std::string, Params> IO::AvailableAttributes(const std::string &variableName,
                                                      const std::string separator,
                                                      const bool fullNameKeys)
{
    helper::CheckForNullptr(m_IO, "in call to IO::AvailableAttributes");
    return m_IO->GetAvailableAttributes(variableName, separator, fullNameKeys);
}

std::string IO::VariableType(const std::string &name) const
{
    helper::CheckForNullptr(m_IO, "in call to IO::VariableType");
    return ToString(m_IO->InquireVariableType(name));
}

std::string IO::AttributeType(const std::string &name) const
{
    helper::CheckForNullptr(m_IO, "in call to IO::AttributeType");
    return ToString(m_IO->InquireAttributeType(name));
}

void IO::AddOperation(const std::string &variable, const std::string &operatorType,
                      const Params &parameters)
{
    helper::CheckForNullptr(m_IO, "in call to IO::AddOperation");
    return m_IO->AddOperation(variable, operatorType, parameters);
}

std::string IO::EngineType() const
{
    helper::CheckForNullptr(m_IO, "in call to IO::EngineType");
    return m_IO->m_EngineType;
}

VariableNT IO::DefineVariable(const DataType type, const std::string &name, const Dims &shape,
                              const Dims &start, const Dims &count, const bool constantDims)
{
    helper::CheckForNullptr(m_IO, "for variable name " + name + ", in call to IO::DefineVariable");
#define declare_type(T)                                                                            \
    if (ToString(type) == GetType<T>())                                                            \
    {                                                                                              \
        return VariableNT(&m_IO->DefineVariable<typename TypeInfo<T>::IOType>(                     \
            name, shape, start, count, constantDims));                                             \
    }
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
    else
    {
        return nullptr;
    }
}

#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
VariableDerived IO::DefineDerivedVariable(const std::string &name, const std::string &expression,
                                          const DerivedVarType varType)
{
    helper::CheckForNullptr(m_IO,
                            "for variable name " + name + ", in call to IO::DefineDerivedVariable");

    return VariableDerived(&m_IO->DefineDerivedVariable(name, expression, varType));
}
#endif
StructDefinition IO::DefineStruct(const std::string &name, const size_t size)
{
    helper::CheckForNullptr(m_IO, "for struct name " + name + ", in call to IO::DefineStruct");
    return StructDefinition(&m_IO->DefineStruct(name, size));
}

VariableNT IO::DefineStructVariable(const std::string &name, const StructDefinition &def,
                                    const Dims &shape, const Dims &start, const Dims &count,
                                    const bool constantDims)
{
    helper::CheckForNullptr(m_IO,
                            "for variable name " + name + ", in call to IO::DefineStructVariable");
    return VariableNT(&m_IO->DefineStructVariable(name, *def.m_StructDefinition, shape, start,
                                                  count, constantDims));
}

VariableNT IO::InquireVariable(const std::string &name)
{
    helper::CheckForNullptr(m_IO, "for variable name " + name + ", in call to IO::InquireVariable");
    auto type = m_IO->InquireVariableType(name);
#define declare_type(T)                                                                            \
    if (type == helper::GetDataType<T>())                                                          \
    {                                                                                              \
        return VariableNT(m_IO->InquireVariable<typename TypeInfo<T>::IOType>(name));              \
    }
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
    else if (type == DataType::Struct)
    {
        return VariableNT(m_IO->InquireStructVariable(name));
    }
    else
    {
        return nullptr;
    }
}

VariableNT IO::InquireStructVariable(const std::string &name)
{
    helper::CheckForNullptr(m_IO,
                            "for variable name " + name + ", in call to IO::InquireStructVariable");

    return VariableNT(m_IO->InquireStructVariable(name));
}

VariableNT IO::InquireStructVariable(const std::string &name, const StructDefinition def)
{
    helper::CheckForNullptr(m_IO,
                            "for variable name " + name + ", in call to IO::InquireStructVariable");

    return VariableNT(m_IO->InquireStructVariable(name, *def.m_StructDefinition));
}

// PRIVATE
IO::IO(core::IO *io) : m_IO(io) {}

// Explicit declaration of the public template methods
// Limits the types
#define declare_template_instantiation(T)                                                          \
    template Variable<T> IO::DefineVariable(const std::string &, const Dims &, const Dims &,       \
                                            const Dims &, const bool);                             \
                                                                                                   \
    template Variable<T> IO::InquireVariable<T>(const std::string &);

ADIOS2_FOREACH_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

#define declare_template_instantiation(T)                                                          \
    template Attribute<T> IO::DefineAttribute(const std::string &, const T *, const size_t,        \
                                              const std::string &, const std::string, const bool); \
                                                                                                   \
    template Attribute<T> IO::DefineAttribute(const std::string &, const T &, const std::string &, \
                                              const std::string, const bool);                      \
                                                                                                   \
    template Attribute<T> IO::InquireAttribute<T>(const std::string &, const std::string &,        \
                                                  const std::string);

ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

std::string ToString(const IO &io) { return std::string("IO(Name: \"" + io.Name() + "\")"); }

} // end namespace adios2
