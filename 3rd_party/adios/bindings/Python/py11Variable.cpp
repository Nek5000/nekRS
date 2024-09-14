/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11Variable.cpp :
 *
 *  Created on: Sep 7, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "py11Variable.h"

#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace py11
{

Variable::Variable(core::VariableBase *variable) : m_VariableBase(variable) {}

Variable::operator bool() const noexcept { return (m_VariableBase == nullptr) ? false : true; }

void Variable::SetShape(const Dims &shape)
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::SetShape");
    m_VariableBase->SetShape(shape);
}

void Variable::SetBlockSelection(const size_t blockID)
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::SetBlockSelection");
    m_VariableBase->SetBlockSelection(blockID);
}

void Variable::SetSelection(const Box<Dims> &selection)
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::SetSelection");
    m_VariableBase->SetSelection(selection);
}

void Variable::SetStepSelection(const Box<size_t> &stepSelection)
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::SetStepSelection");
    m_VariableBase->SetStepSelection(stepSelection);
}

size_t Variable::SelectionSize() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::SelectionSize");

    const adios2::DataType typeCpp = m_VariableBase->m_Type;
    size_t size = 0;

    if (typeCpp == adios2::DataType::Struct)
    {
        // not supported
    }
#define declare_template_instantiation(T)                                                          \
    else if (typeCpp == adios2::helper::GetDataType<T>())                                          \
    {                                                                                              \
        const adios2::core::Variable<T> *variable =                                                \
            dynamic_cast<const adios2::core::Variable<T> *>(m_VariableBase);                       \
        size = variable->SelectionSize();                                                          \
    }
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

    return size;
}

size_t Variable::AddOperation(const Operator op, const Params &parameters)
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::AddOperation");
    auto params = op.Parameters();
    for (const auto &p : parameters)
    {
        params[p.first] = p.second;
    }
    return m_VariableBase->AddOperation(op.m_Type, params);
}

size_t Variable::AddOperation(const std::string &op, const Params &parameters)
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::AddOperation");
    return m_VariableBase->AddOperation(op, parameters);
}

std::vector<Operator> Variable::Operations() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::Operations");
    std::vector<Operator> operations;
    operations.reserve(m_VariableBase->m_Operations.size());

    for (const auto &op : m_VariableBase->m_Operations)
    {
        operations.push_back(Operator(op->m_TypeString, &op->GetParameters()));
    }
    return operations;
}

void Variable::RemoveOperations() { m_VariableBase->RemoveOperations(); }

std::string Variable::Name() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::Name");
    return m_VariableBase->m_Name;
}

std::string Variable::Type() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::Type");
    return ToString(m_VariableBase->m_Type);
}

size_t Variable::Sizeof() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::Sizeof");
    return m_VariableBase->m_ElementSize;
}

adios2::ShapeID Variable::ShapeID() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::ShapeID");
    return m_VariableBase->m_ShapeID;
}

Dims Variable::Shape(const size_t step) const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::Shape");

    const adios2::DataType typeCpp = m_VariableBase->m_Type;
    Dims shape;

    if (typeCpp == adios2::DataType::Struct)
    {
        // not supported
    }
#define declare_template_instantiation(T)                                                          \
    else if (typeCpp == adios2::helper::GetDataType<T>())                                          \
    {                                                                                              \
        const adios2::core::Variable<T> *variable =                                                \
            dynamic_cast<const adios2::core::Variable<T> *>(m_VariableBase);                       \
        shape = variable->Shape(step);                                                             \
    }
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

    return shape;
}

Dims Variable::Start() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::Start");
    return m_VariableBase->m_Start;
}

Dims Variable::Count() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::Count");

    const adios2::DataType typeCpp = m_VariableBase->m_Type;
    Dims count;

    if (typeCpp == adios2::DataType::Struct)
    {
        // not supported
    }
#define declare_template_instantiation(T)                                                          \
    else if (typeCpp == adios2::helper::GetDataType<T>())                                          \
    {                                                                                              \
        const adios2::core::Variable<T> *variable =                                                \
            dynamic_cast<const adios2::core::Variable<T> *>(m_VariableBase);                       \
        count = variable->Count();                                                                 \
    }
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

    return count;
}

size_t Variable::Steps() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::Steps");
    return m_VariableBase->m_AvailableStepsCount;
}

size_t Variable::StepsStart() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::StepsStart");
    return m_VariableBase->m_AvailableStepsStart;
}

size_t Variable::BlockID() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::BlockID");
    return m_VariableBase->m_BlockID;
}

size_t Variable::SingleValue() const
{
    helper::CheckForNullptr(m_VariableBase, "in call to Variable::SingleValue");
    return m_VariableBase->m_SingleValue;
}

// size_t Variable::AddOperation(const Operator op, const Params &parameters) {}

// std::vector<Operation> Variable::Operations() const {}

} // end namespace py11
} // end namespace adios2
