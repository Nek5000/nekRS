/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Stream.tcc
 *
 *  Created on: Jan 5, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_STREAM_TCC_
#define ADIOS2_CORE_STREAM_TCC_

#include "Stream.h"

#include "adios2/core/Variable.h"
#include "adios2/helper/adiosLog.h"

namespace adios2
{
namespace core
{

template <class T>
void Stream::WriteAttribute(const std::string &name, const T &value,
                            const std::string &variableName, const std::string separator,
                            const bool endStep)
{
    m_IO->DefineAttribute<T>(name, value, variableName, separator);
    CheckOpen();
    if (!m_StepStatus)
    {
        m_Engine->BeginStep();
        m_StepStatus = true;
    }

    if (endStep)
    {
        m_Engine->EndStep();
        m_StepStatus = false;
    }
}

template <class T>
void Stream::WriteAttribute(const std::string &name, const T *array, const size_t elements,
                            const std::string &variableName, const std::string separator,
                            const bool endStep)
{
    m_IO->DefineAttribute<T>(name, array, elements, variableName, separator);
    CheckOpen();
    if (!m_StepStatus)
    {
        m_Engine->BeginStep();
        m_StepStatus = true;
    }

    if (endStep)
    {
        m_Engine->EndStep();
        m_StepStatus = false;
    }
}

template <class T>
void Stream::Write(const std::string &name, const T *data, const Dims &shape, const Dims &start,
                   const Dims &count, const adios2::vParams &operations, const bool endStep)
{
    Variable<T> *variable = m_IO->InquireVariable<T>(name);

    if (variable == nullptr)
    {
        variable = &m_IO->DefineVariable<T>(name, shape, start, count, false);
    }
    else
    {
        if (!shape.empty() && !variable->m_SingleValue)
        {
            variable->SetShape(shape);
        }

        if (!start.empty() && !count.empty())
        {
            variable->SetSelection(Box<Dims>(start, count));
        }
    }

    CheckOpen();
    if (!m_StepStatus)
    {
        m_Engine->BeginStep();
        m_StepStatus = true;
    }

    // boiler plate for operations
    if (!operations.empty())
    {
        variable->m_Operations.clear();
        for (const auto &operation : operations)
        {
            variable->AddOperation(operation.first, operation.second);
        }
    }

    m_Engine->Put(*variable, data, adios2::Mode::Sync);

    if (endStep)
    {
        m_Engine->EndStep();
        m_StepStatus = false;
    }
}

template <class T>
void Stream::Write(const std::string &name, const T &datum, const bool isLocal, const bool endStep)
{
    const T datumLocal = datum;
    if (isLocal)
    {
        Write(name, &datumLocal, {static_cast<size_t>(adios2::LocalValueDim)}, {}, {}, vParams(),
              endStep);
    }
    else
    {
        Write(name, &datumLocal, {}, {}, {}, vParams(), endStep);
    }
}

template <class T>
void Stream::Read(const std::string &name, T *values, const size_t blockID)
{
    CheckPCommon(name, values);

    Variable<T> *variable = m_IO->InquireVariable<T>(name);
    if (variable == nullptr)
    {
        values = nullptr;
        return;
    }

    SetBlockSelectionCommon(*variable, blockID);
    GetPCommon(*variable, values);
}

template <class T>
void Stream::Read(const std::string &name, T *values, const Box<size_t> &stepSelection,
                  const size_t blockID)
{
    CheckPCommon(name, values);

    Variable<T> *variable = m_IO->InquireVariable<T>(name);
    if (variable == nullptr)
    {
        values = nullptr;
        return;
    }

    SetBlockSelectionCommon(*variable, blockID);
    variable->SetStepSelection(stepSelection);
    GetPCommon(*variable, values);
}

template <class T>
void Stream::Read(const std::string &name, T *values, const Box<Dims> &selection,
                  const size_t blockID)
{
    CheckPCommon(name, values);

    Variable<T> *variable = m_IO->InquireVariable<T>(name);
    if (variable == nullptr)
    {
        values = nullptr;
        return;
    }

    SetBlockSelectionCommon(*variable, blockID);
    variable->SetSelection(selection);
    GetPCommon(*variable, values);
}

template <class T>
void Stream::Read(const std::string &name, T *values, const Box<Dims> &selection,
                  const Box<size_t> &stepSelection, const size_t blockID)
{
    CheckPCommon(name, values);

    Variable<T> *variable = m_IO->InquireVariable<T>(name);
    if (variable == nullptr)
    {
        values = nullptr;
        return;
    }

    SetBlockSelectionCommon(*variable, blockID);
    variable->SetSelection(selection);
    variable->SetStepSelection(stepSelection);
    GetPCommon(*variable, values);
}

template <class T>
std::vector<T> Stream::Read(const std::string &name, const size_t blockID)
{
    Variable<T> *variable = m_IO->InquireVariable<T>(name);
    if (variable == nullptr)
    {
        return std::vector<T>();
    }
    SetBlockSelectionCommon(*variable, blockID);
    return GetCommon(*variable);
}

template <class T>
std::vector<T> Stream::Read(const std::string &name, const Box<size_t> &stepsSelection,
                            const size_t blockID)
{
    Variable<T> *variable = m_IO->InquireVariable<T>(name);
    if (variable == nullptr)
    {
        return std::vector<T>();
    }
    SetBlockSelectionCommon(*variable, blockID);
    variable->SetStepSelection(stepsSelection);
    return GetCommon(*variable);
}

template <class T>
std::vector<T> Stream::Read(const std::string &name, const Box<Dims> &selection,
                            const size_t blockID)
{
    Variable<T> *variable = m_IO->InquireVariable<T>(name);
    if (variable == nullptr)
    {
        return std::vector<T>();
    }

    SetBlockSelectionCommon(*variable, blockID);
    variable->SetSelection(selection);
    return GetCommon(*variable);
}

template <class T>
std::vector<T> Stream::Read(const std::string &name, const Box<Dims> &selection,
                            const Box<size_t> &stepSelection, const size_t blockID)
{
    Variable<T> *variable = m_IO->InquireVariable<T>(name);
    if (variable == nullptr)
    {
        return std::vector<T>();
    }

    SetBlockSelectionCommon(*variable, blockID);
    variable->SetSelection(selection);
    variable->SetStepSelection(stepSelection);
    return GetCommon(*variable);
}

template <class T>
void Stream::ReadAttribute(const std::string &name, T *data, const std::string &variableName,
                           const std::string separator)
{
    Attribute<T> *attribute = m_IO->InquireAttribute<T>(name, variableName, separator);

    if (attribute == nullptr)
    {
        return;
    }

    if (attribute->m_IsSingleValue)
    {
        data[0] = attribute->m_DataSingleValue;
    }
    else
    {
        std::copy(attribute->m_DataArray.begin(), attribute->m_DataArray.end(), data);
    }
}

// PRIVATE
template <class T>
std::vector<T> Stream::GetCommon(Variable<T> &variable)
{
    try
    {
        std::vector<T> values(variable.SelectionSize());
        CheckOpen();
        m_Engine->Get(variable, values.data(), adios2::Mode::Sync);
        return values;
    }
    catch (std::exception &e)
    {
        helper::ThrowNested<std::runtime_error>("Core", "Stream", "GetCommon",
                                                "couldn't read variable " + variable.m_Name + ": " +
                                                    e.what());
    }
    return std::vector<T>();
}

template <class T>
void Stream::GetPCommon(Variable<T> &variable, T *values)
{
    try
    {
        CheckOpen();
        m_Engine->Get(variable, values, adios2::Mode::Sync);
    }
    catch (std::exception &e)
    {
        helper::ThrowNested<std::runtime_error>("Core", "Stream", "GetCommon",
                                                "couldn't read pointer variable " +
                                                    variable.m_Name + ": " + e.what());
    }
}

template <class T>
void Stream::CheckPCommon(const std::string &name, const T *values) const
{
    if (values == nullptr)
    {
        helper::Throw<std::runtime_error>("Core", "Stream", "CheckPCommon",
                                          "passed null values pointer for variable " + name +
                                              ", in call to read pointer");
    }
}

template <class T>
void Stream::SetBlockSelectionCommon(Variable<T> &variable, const size_t blockID)
{
    if (variable.m_ShapeID == ShapeID::LocalArray)
    {
        variable.SetBlockSelection(blockID);
    }
    else
    {
        if (blockID != 0)
        {
            helper::Throw<std::invalid_argument>("Core", "Stream", "SetBlockSelectionCommon",
                                                 "in variable " + variable.m_Name +
                                                     " only set blockID > 0 for variables "
                                                     "with ShapeID::LocalArray, in call to read");
        }
    }
}

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_STREAM_TCC_ */
