/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Variable.cpp :
 *
 *  Created on: Apr 18, 2022
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "VariableNT.h"
#include "Types.h"
#include "adios2/core/VariableBase.h"
#include "adios2/helper/adiosFunctions.h"

namespace adios2
{

StructDefinition::StructDefinition(core::StructDefinition *ptr) : m_StructDefinition(ptr) {}

void StructDefinition::AddField(const std::string &name, const size_t offset, const DataType type,
                                const size_t size)
{
    helper::CheckForNullptr(m_StructDefinition, "in call to StructDefinition::AddField");
    m_StructDefinition->AddField(name, offset, type, size);
}

void StructDefinition::Freeze() noexcept
{
    helper::CheckForNullptr(m_StructDefinition, "in call to StructDefinition::Freeze");
    m_StructDefinition->Freeze();
}

size_t StructDefinition::StructSize() const noexcept
{
    helper::CheckForNullptr(m_StructDefinition, "in call to StructDefinition::StructSize");
    return m_StructDefinition->StructSize();
}

std::string StructDefinition::StructName() const noexcept
{
    helper::CheckForNullptr(m_StructDefinition, "in call to StructDefinition::StructName");
    return m_StructDefinition->StructName();
}

size_t StructDefinition::Fields() const noexcept
{
    helper::CheckForNullptr(m_StructDefinition, "in call to StructDefinition::Fields");
    return m_StructDefinition->Fields();
}
std::string StructDefinition::Name(const size_t index) const
{
    helper::CheckForNullptr(m_StructDefinition, "in call to StructDefinition::Name");
    return m_StructDefinition->Name(index);
}
size_t StructDefinition::Offset(const size_t index) const
{
    helper::CheckForNullptr(m_StructDefinition, "in call to StructDefinition::Offset");
    return m_StructDefinition->Offset(index);
}
DataType StructDefinition::Type(const size_t index) const
{
    helper::CheckForNullptr(m_StructDefinition, "in call to StructDefinition::Type");
    return m_StructDefinition->Type(index);
}
size_t StructDefinition::ElementCount(const size_t index) const
{
    helper::CheckForNullptr(m_StructDefinition, "in call to StructDefinition::ElementCount");
    return m_StructDefinition->ElementCount(index);
}

VariableNT::VariableNT(core::VariableBase *variable) : m_Variable(variable) {}

VariableNT::operator bool() const noexcept { return m_Variable != nullptr; }

void VariableNT::SetMemorySpace(const MemorySpace mem)
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::SetMemorySpace");
    m_Variable->SetMemorySpace(mem);
}

void VariableNT::SetShape(const Dims &shape)
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::SetShape");
    m_Variable->SetShape(shape);
}

void VariableNT::SetBlockSelection(const size_t blockID)
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::SetBlockSelection");
    m_Variable->SetBlockSelection(blockID);
}

void VariableNT::SetSelection(const Box<Dims> &selection)
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::SetSelection");
    m_Variable->SetSelection(selection);
}

void VariableNT::SetMemorySelection(const Box<Dims> &memorySelection)
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::SetMemorySelection");
    m_Variable->SetMemorySelection(memorySelection);
}

void VariableNT::SetStepSelection(const Box<size_t> &stepSelection)
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::SetStepSelection");
    m_Variable->SetStepSelection(stepSelection);
}

size_t VariableNT::SelectionSize() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::SelectionSize");
    auto type = ToString(m_Variable->m_Type);
#define declare_type(T)                                                                            \
    if (type == GetType<T>())                                                                      \
    {                                                                                              \
        return reinterpret_cast<core::Variable<T> *>(m_Variable)->SelectionSize();                 \
    }
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
    helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "SelectionSize",
                                      "invalid data type " + type);
    return 0;
}

std::string VariableNT::Name() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Name");
    return m_Variable->m_Name;
}

std::string VariableNT::Type() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Type");
    return ToString(m_Variable->m_Type);
}

size_t VariableNT::Sizeof() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Sizeof");
    return m_Variable->m_ElementSize;
}

adios2::ShapeID VariableNT::ShapeID() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::ShapeID");
    return m_Variable->m_ShapeID;
}

Dims VariableNT::Shape(const size_t step) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Shape");
    return m_Variable->Shape(step);
}

Dims VariableNT::Start() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Start");
    return m_Variable->m_Start;
}

Dims VariableNT::Count() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Count");
    auto type = m_Variable->m_Type;
#define declare_type(T)                                                                            \
    if (type == helper::GetDataType<T>())                                                          \
    {                                                                                              \
        return reinterpret_cast<core::Variable<T> *>(m_Variable)->Count();                         \
    }
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
    else if (type == DataType::Struct)
    {
        return reinterpret_cast<core::VariableStruct *>(m_Variable)->m_Count;
    }
    helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "Count",
                                      "invalid data type " + ToString(type));
    return Dims();
}

size_t VariableNT::Steps() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Steps");
    return m_Variable->m_AvailableStepsCount;
}

size_t VariableNT::StepsStart() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::StepsStart");
    return m_Variable->m_AvailableStepsStart;
}

size_t VariableNT::BlockID() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::BlockID");
    return m_Variable->m_BlockID;
}

size_t VariableNT::AddOperation(const Operator op, const Params &parameters)
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::AddOperation");
    if (!op)
    {
        helper::Throw<std::invalid_argument>("bindings::CXX11", "VariableNT", "AddOperation",
                                             "invalid operation");
    }
    auto params = op.Parameters();
    for (const auto &p : parameters)
    {
        params[p.first] = p.second;
    }
    return m_Variable->AddOperation(op.m_Type, params);
}

size_t VariableNT::AddOperation(const std::string &type, const Params &parameters)
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::AddOperation");
    return m_Variable->AddOperation(type, parameters);
}

std::vector<Operator> VariableNT::Operations() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Operations");
    std::vector<Operator> operations;
    operations.reserve(m_Variable->m_Operations.size());
    for (const auto &op : m_Variable->m_Operations)
    {
        operations.push_back(Operator(op->m_TypeString, &op->GetParameters()));
    }
    return operations;
}

void VariableNT::RemoveOperations()
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::RemoveOperations");
    m_Variable->RemoveOperations();
}

size_t VariableNT::StructFields() const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::StructFields");
    if (m_Variable->m_Type != DataType::Struct)
    {
        helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "StructFields",
                                          "invalid data type " + ToString(m_Variable->m_Type) +
                                              ", only Struct type supports this API");
    }
    if (reinterpret_cast<core::VariableStruct *>(m_Variable)->m_ReadStructDefinition)
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_ReadStructDefinition->Fields();
    else
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_WriteStructDefinition->Fields();
}
std::string VariableNT::StructFieldName(const size_t index) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::StructFieldName");
    if (m_Variable->m_Type != DataType::Struct)
    {
        helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "StructFieldName",
                                          "invalid data type " + ToString(m_Variable->m_Type) +
                                              ", only Struct type supports this API");
    }
    if (reinterpret_cast<core::VariableStruct *>(m_Variable)->m_ReadStructDefinition)
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_ReadStructDefinition->Name(index);
    else
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_WriteStructDefinition->Name(index);
}
size_t VariableNT::StructFieldOffset(const size_t index) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::StructFieldOffset");
    if (m_Variable->m_Type != DataType::Struct)
    {
        helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "StructFieldOffset",
                                          "invalid data type " + ToString(m_Variable->m_Type) +
                                              ", only Struct type supports this API");
    }
    if (reinterpret_cast<core::VariableStruct *>(m_Variable)->m_ReadStructDefinition)
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_ReadStructDefinition->Offset(index);
    else
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_WriteStructDefinition->Offset(index);
}
DataType VariableNT::StructFieldType(const size_t index) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::StructFieldType");
    if (m_Variable->m_Type != DataType::Struct)
    {
        helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "StructFieldType",
                                          "invalid data type " + ToString(m_Variable->m_Type) +
                                              ", only Struct type supports this API");
    }
    if (reinterpret_cast<core::VariableStruct *>(m_Variable)->m_ReadStructDefinition)
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_ReadStructDefinition->Type(index);
    else
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_WriteStructDefinition->Type(index);
}
size_t VariableNT::StructFieldElementCount(const size_t index) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::StructFieldElementCount");
    if (m_Variable->m_Type != DataType::Struct)
    {
        helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT",
                                          "StructFieldElementCount",
                                          "invalid data type " + ToString(m_Variable->m_Type) +
                                              ", only Struct type supports this API");
    }
    if (reinterpret_cast<core::VariableStruct *>(m_Variable)->m_ReadStructDefinition)
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_ReadStructDefinition->ElementCount(index);
    else
        return reinterpret_cast<core::VariableStruct *>(m_Variable)
            ->m_WriteStructDefinition->ElementCount(index);
}

std::pair<VariableNT::T, VariableNT::T> VariableNT::MinMax(const size_t step) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::MinMax");
    return {Min(step), Max(step)};
}

VariableNT::T VariableNT::Min(const size_t step) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Min");
    if (m_Variable->m_Type == DataType::Int8)
    {
        VariableNT::T ret = {0};
        ret.Int8 = reinterpret_cast<core::Variable<int8_t> *>(m_Variable)->Min(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::UInt8)
    {
        VariableNT::T ret = {0};
        ret.UInt8 = reinterpret_cast<core::Variable<uint8_t> *>(m_Variable)->Min(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Int16)
    {
        VariableNT::T ret = {0};
        ret.Int16 = reinterpret_cast<core::Variable<int16_t> *>(m_Variable)->Min(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::UInt16)
    {
        VariableNT::T ret = {0};
        ret.UInt16 = reinterpret_cast<core::Variable<uint16_t> *>(m_Variable)->Min(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Int32)
    {
        VariableNT::T ret = {0};
        ret.Int32 = reinterpret_cast<core::Variable<int32_t> *>(m_Variable)->Min(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::UInt32)
    {
        VariableNT::T ret = {0};
        ret.UInt32 = reinterpret_cast<core::Variable<uint32_t> *>(m_Variable)->Min(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Int64)
    {
        VariableNT::T ret = {0};
        ret.Int64 = reinterpret_cast<core::Variable<int64_t> *>(m_Variable)->Min(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::UInt64)
    {
        VariableNT::T ret = {0};
        ret.UInt64 = reinterpret_cast<core::Variable<uint64_t> *>(m_Variable)->Min(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Float)
    {
        VariableNT::T ret = {0};
        ret.Float = reinterpret_cast<core::Variable<float> *>(m_Variable)->Min(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Double)
    {
        VariableNT::T ret = {0};
        ret.Double = reinterpret_cast<core::Variable<double> *>(m_Variable)->Min(step);
        return ret;
    }
    helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "Min",
                                      "invalid data type " + ToString(m_Variable->m_Type) +
                                          ", only basic numeric types support this API");
    return {0};
}

VariableNT::T VariableNT::Max(const size_t step) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::Max");
    if (m_Variable->m_Type == DataType::Int8)
    {
        VariableNT::T ret = {0};
        ret.Int8 = reinterpret_cast<core::Variable<int8_t> *>(m_Variable)->Max(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::UInt8)
    {
        VariableNT::T ret = {0};
        ret.UInt8 = reinterpret_cast<core::Variable<uint8_t> *>(m_Variable)->Max(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Int16)
    {
        VariableNT::T ret = {0};
        ret.Int16 = reinterpret_cast<core::Variable<int16_t> *>(m_Variable)->Max(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::UInt16)
    {
        VariableNT::T ret = {0};
        ret.UInt16 = reinterpret_cast<core::Variable<uint16_t> *>(m_Variable)->Max(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Int32)
    {
        VariableNT::T ret = {0};
        ret.Int32 = reinterpret_cast<core::Variable<int32_t> *>(m_Variable)->Max(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::UInt32)
    {
        VariableNT::T ret = {0};
        ret.UInt32 = reinterpret_cast<core::Variable<uint32_t> *>(m_Variable)->Max(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Int64)
    {
        VariableNT::T ret = {0};
        ret.Int64 = reinterpret_cast<core::Variable<int64_t> *>(m_Variable)->Max(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::UInt64)
    {
        VariableNT::T ret = {0};
        ret.UInt64 = reinterpret_cast<core::Variable<uint64_t> *>(m_Variable)->Max(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Float)
    {
        VariableNT::T ret = {0};
        ret.Float = reinterpret_cast<core::Variable<float> *>(m_Variable)->Max(step);
        return ret;
    }
    else if (m_Variable->m_Type == DataType::Double)
    {
        VariableNT::T ret = {0};
        ret.Double = reinterpret_cast<core::Variable<double> *>(m_Variable)->Max(step);
        return ret;
    }
    helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "Max",
                                      "invalid data type " + ToString(m_Variable->m_Type) +
                                          ", only basic numeric types support this API");
    return {0};
}

double VariableNT::MinDouble(const size_t step) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::MinDouble");
    if (m_Variable->m_Type == DataType::Int8)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<int8_t> *>(m_Variable)->Min(step));
    }
    else if (m_Variable->m_Type == DataType::UInt8)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<uint8_t> *>(m_Variable)->Min(step));
    }
    else if (m_Variable->m_Type == DataType::Int16)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<int16_t> *>(m_Variable)->Min(step));
    }
    else if (m_Variable->m_Type == DataType::UInt16)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<uint16_t> *>(m_Variable)->Min(step));
    }
    else if (m_Variable->m_Type == DataType::Int32)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<int32_t> *>(m_Variable)->Min(step));
    }
    else if (m_Variable->m_Type == DataType::UInt32)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<uint32_t> *>(m_Variable)->Min(step));
    }
    else if (m_Variable->m_Type == DataType::Int64)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<int64_t> *>(m_Variable)->Min(step));
    }
    else if (m_Variable->m_Type == DataType::UInt64)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<uint64_t> *>(m_Variable)->Min(step));
    }
    else if (m_Variable->m_Type == DataType::Float)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<float> *>(m_Variable)->Min(step));
    }
    else if (m_Variable->m_Type == DataType::Double)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<double> *>(m_Variable)->Min(step));
    }
    helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "MinDouble",
                                      "invalid data type " + ToString(m_Variable->m_Type) +
                                          ", only basic numeric types support this API");
    return 0;
}

double VariableNT::MaxDouble(const size_t step) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::MaxDouble");
    if (m_Variable->m_Type == DataType::Int8)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<int8_t> *>(m_Variable)->Max(step));
    }
    else if (m_Variable->m_Type == DataType::UInt8)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<uint8_t> *>(m_Variable)->Max(step));
    }
    else if (m_Variable->m_Type == DataType::Int16)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<int16_t> *>(m_Variable)->Max(step));
    }
    else if (m_Variable->m_Type == DataType::UInt16)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<uint16_t> *>(m_Variable)->Max(step));
    }
    else if (m_Variable->m_Type == DataType::Int32)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<int32_t> *>(m_Variable)->Max(step));
    }
    else if (m_Variable->m_Type == DataType::UInt32)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<uint32_t> *>(m_Variable)->Max(step));
    }
    else if (m_Variable->m_Type == DataType::Int64)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<int64_t> *>(m_Variable)->Max(step));
    }
    else if (m_Variable->m_Type == DataType::UInt64)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<uint64_t> *>(m_Variable)->Max(step));
    }
    else if (m_Variable->m_Type == DataType::Float)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<float> *>(m_Variable)->Max(step));
    }
    else if (m_Variable->m_Type == DataType::Double)
    {
        return static_cast<double>(
            reinterpret_cast<core::Variable<double> *>(m_Variable)->Max(step));
    }
    helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT", "MaxDouble",
                                      "invalid data type " + ToString(m_Variable->m_Type) +
                                          ", only basic numeric types support this API");
    return 0;
}

std::pair<double, double> VariableNT::MinMaxDouble(const size_t step) const
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::MinMaxDouble");
    return {MinDouble(step), MaxDouble(step)};
}

StructDefinition VariableNT::GetWriteStructDef() noexcept
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::StructFieldElementCount");
    if (m_Variable->m_Type != DataType::Struct)
    {
        helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT",
                                          "StructFieldElementCount",
                                          "invalid data type " + ToString(m_Variable->m_Type) +
                                              ", only Struct type supports this API");
    }
    core::StructDefinition *CoreSD =
        reinterpret_cast<core::VariableStruct *>(m_Variable)->GetWriteStructDef();
    return StructDefinition(CoreSD);
}

StructDefinition VariableNT::GetReadStructDef() noexcept
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::StructFieldElementCount");
    if (m_Variable->m_Type != DataType::Struct)
    {
        helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT",
                                          "StructFieldElementCount",
                                          "invalid data type " + ToString(m_Variable->m_Type) +
                                              ", only Struct type supports this API");
    }
    auto CoreSD = reinterpret_cast<core::VariableStruct *>(m_Variable)->GetReadStructDef();
    return StructDefinition(CoreSD);
}

void VariableNT::SetReadStructDef(const StructDefinition &def)
{
    helper::CheckForNullptr(m_Variable, "in call to VariableNT::StructFieldElementCount");
    if (m_Variable->m_Type != DataType::Struct)
    {
        helper::Throw<std::runtime_error>("bindings::CXX11", "VariableNT",
                                          "StructFieldElementCount",
                                          "invalid data type " + ToString(m_Variable->m_Type) +
                                              ", only Struct type supports this API");
    }
    reinterpret_cast<core::VariableStruct *>(m_Variable)->SetReadStructDef(def.m_StructDefinition);
}

} // end namespace adios2
