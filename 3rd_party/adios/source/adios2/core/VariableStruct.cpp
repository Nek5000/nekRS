/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * VariableStruct.cpp
 *
 *  Created on: Apr 24, 2022
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "VariableStruct.h"
#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosType.h"

namespace adios2
{
namespace core
{

StructDefinition::StructDefinition(const std::string &name, const size_t size)
: m_Name(name), m_StructSize(size)
{
}

void StructDefinition::AddField(const std::string &name, const size_t offset, const DataType type,
                                const size_t elementcount)
{
    if (m_Frozen)
    {
        helper::Throw<std::runtime_error>("core", "VariableStruct::StructDefinition", "AddField",
                                          "struct definition already frozen");
    }
    if (type == DataType::None || type == DataType::Struct)
    {
        helper::Throw<std::invalid_argument>("core", "VariableStruct::StructDefinition", "AddField",
                                             "invalid data type");
    }
    if (offset + helper::GetDataTypeSize(type) * elementcount > m_StructSize)
    {
        helper::Throw<std::runtime_error>("core", "VariableStruct::StructDefinition", "AddField",
                                          "exceeded struct size");
    }
    m_Definition.emplace_back();
    auto &d = m_Definition.back();
    d.Name = name;
    d.Offset = offset;
    d.Type = type;
    d.ElementCount = elementcount;
}

void StructDefinition::Freeze() noexcept { m_Frozen = true; }

size_t StructDefinition::StructSize() const noexcept { return m_StructSize; }

size_t StructDefinition::Fields() const noexcept { return m_Definition.size(); }

std::string StructDefinition::StructName() const noexcept { return m_Name; }

std::string StructDefinition::Name(const size_t index) const
{
    if (index >= m_Definition.size())
    {
        helper::Throw<std::invalid_argument>("core", "VariableStruct::StructDefinition", "Name",
                                             "invalid index");
    }
    return m_Definition[index].Name;
}

size_t StructDefinition::Offset(const size_t index) const
{
    if (index >= m_Definition.size())
    {
        helper::Throw<std::invalid_argument>("core", "VariableStruct::StructDefinition", "Offset",
                                             "invalid index");
    }
    return m_Definition[index].Offset;
}

DataType StructDefinition::Type(const size_t index) const
{
    if (index >= m_Definition.size())
    {
        helper::Throw<std::invalid_argument>("core", "VariableStruct::StructDefinition", "Type",
                                             "invalid index");
    }
    return m_Definition[index].Type;
}

size_t StructDefinition::ElementCount(const size_t index) const
{
    if (index >= m_Definition.size())
    {
        helper::Throw<std::invalid_argument>("core", "VariableStruct::StructDefinition",
                                             "ElementCount", "invalid index");
    }
    return m_Definition[index].ElementCount;
}

VariableStruct::VariableStruct(const std::string &name, const StructDefinition &def,
                               const Dims &shape, const Dims &start, const Dims &count,
                               const bool constantDims)
: VariableBase(name, DataType::Struct, def.StructSize(), shape, start, count, constantDims)
{
    m_WriteStructDefinition = const_cast<StructDefinition *>(&def);
    m_ReadStructDefinition = nullptr;
}

void VariableStruct::SetData(const void *data) noexcept { m_Data = const_cast<void *>(data); }

void *VariableStruct::GetData() const noexcept { return m_Data; }

StructDefinition *VariableStruct::GetWriteStructDef() noexcept { return m_WriteStructDefinition; }

StructDefinition *VariableStruct::GetReadStructDef() noexcept { return m_ReadStructDefinition; }

void VariableStruct::SetReadStructDef(const StructDefinition *def)
{
    m_ReadStructDefinition = (StructDefinition *)def;
}

} // end namespace core
} // end namespace adios2
