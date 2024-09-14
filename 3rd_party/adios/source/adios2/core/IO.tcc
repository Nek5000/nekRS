/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO.tcc template implementations with fix types and specializations
 *
 *  Created on: May 15, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_IO_TCC_
#define ADIOS2_CORE_IO_TCC_

#include "IO.h"

/// \cond EXCLUDE_FROM_DOXYGEN
#include <iostream>
#include <memory>
#include <stdexcept> //std::invalid_argument
/// \endcond

#include "adios2/common/ADIOSMacros.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/helper/adiosType.h"
#include <adios2-perfstubs-interface.h>
#include <adios2/core/Engine.h>

namespace adios2
{
namespace core
{

template <class T>
Variable<T> &IO::DefineVariable(const std::string &name, const Dims &shape, const Dims &start,
                                const Dims &count, const bool constantDims)
{
    PERFSTUBS_SCOPED_TIMER("IO::DefineVariable");

    {
        auto itVariable = m_Variables.find(name);
        if (itVariable != m_Variables.end())
        {
            helper::Throw<std::invalid_argument>("Core", "IO", "DefineVariable",
                                                 "variable " + name + " already defined in IO " +
                                                     m_Name);
        }
    }

    auto itVariablePair = m_Variables.emplace(name, std::unique_ptr<VariableBase>(new Variable<T>(
                                                        name, shape, start, count, constantDims)));

    Variable<T> &variable = static_cast<Variable<T> &>(*itVariablePair.first->second);

    // check IO placeholder for variable operations
    auto itOperations = m_VarOpsPlaceholder.find(name);
    if (itOperations != m_VarOpsPlaceholder.end())
    {
        variable.m_Operations.reserve(itOperations->second.size());
        for (auto &operation : itOperations->second)
        {
            variable.AddOperation(operation.first, operation.second);
        }
    }

#if defined(ADIOS2_HAVE_KOKKOS) || defined(ADIOS2_HAVE_GPU_SUPPORT)
    variable.m_BaseLayout = m_ArrayOrder;
#endif
    return variable;
}

template <class T>
Variable<T> *IO::InquireVariable(const std::string &name) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::InquireVariable");
    auto itVariable = m_Variables.find(name);

    if (m_Variables.empty())
    {
        for (auto &e : m_Engines)
        {
            e.second->NotifyEngineNoVarsQuery();
        }
    }
    if (itVariable == m_Variables.end())
    {
        return nullptr;
    }

    if (itVariable->second->m_Type != helper::GetDataType<T>())
    {
        return nullptr;
    }

    Variable<T> *variable = static_cast<Variable<T> *>(itVariable->second.get());
    if (m_ReadStreaming)
    {
        if (!variable->IsValidStep(m_EngineStep + 1))
        {
            return nullptr;
        }
    }
    return variable;
}

template <class T>
Attribute<T> &IO::DefineAttribute(const std::string &name, const T &value,
                                  const std::string &variableName, const std::string separator,
                                  const bool allowModification)
{
    PERFSTUBS_SCOPED_TIMER("IO::DefineAttribute");
    if (!variableName.empty() && InquireVariableType(variableName) == DataType::None)
    {
        helper::Throw<std::invalid_argument>("Core", "IO", "DefineAttribute",
                                             "variable " + variableName +
                                                 " doesn't exist, can't associate attribute " +
                                                 name + ", in call to DefineAttribute");
    }

    const std::string globalName = helper::GlobalName(name, variableName, separator);

    auto itExistingAttribute = m_Attributes.find(globalName);
    if (itExistingAttribute != m_Attributes.end())
    {
        if (itExistingAttribute->second->m_Type == helper::GetDataType<T>())
        {
            if (!itExistingAttribute->second->Equals(static_cast<const void *>(&value), 1))
            {
                Attribute<T> &a = static_cast<Attribute<T> &>(*itExistingAttribute->second);

                a.Modify(value);
                void *Data = &a.m_DataSingleValue;
                if (a.m_DataArray.size() != 0)
                    Data = a.m_DataArray.data();
                for (auto &e : m_Engines)
                {
                    e.second->NotifyEngineAttribute(globalName, itExistingAttribute->second.get(),
                                                    Data);
                }
            }
        }
        else
        {
            helper::Throw<std::invalid_argument>(
                "Core", "IO", "DefineAttribute",
                "modifiable attribute " + globalName + " has been defined with type " +
                    ToString(itExistingAttribute->second->m_Type) + ". Type cannot be changed to " +
                    ToString(helper::GetDataType<T>()));
        }
        return static_cast<Attribute<T> &>(*itExistingAttribute->second);
    }
    else
    {
        auto itAttributePair = m_Attributes.emplace(
            globalName,
            std::unique_ptr<AttributeBase>(new Attribute<T>(globalName, value, allowModification)));
        for (auto &e : m_Engines)
        {
            Attribute<T> &a = static_cast<Attribute<T> &>(*itAttributePair.first->second);
            void *Data = &a.m_DataSingleValue;
            if (a.m_DataArray.size() != 0)
                Data = a.m_DataArray.data();

            e.second->NotifyEngineAttribute(globalName, itAttributePair.first->second.get(), Data);
        }
        return static_cast<Attribute<T> &>(*itAttributePair.first->second);
    }
}

template <class T>
Attribute<T> &IO::DefineAttribute(const std::string &name, const T *array, const size_t elements,
                                  const std::string &variableName, const std::string separator,
                                  const bool allowModification)
{
    PERFSTUBS_SCOPED_TIMER("IO::DefineAttribute");
    if (!variableName.empty() && InquireVariableType(variableName) == DataType::None)
    {
        helper::Throw<std::invalid_argument>("Core", "IO", "DefineAttribute",
                                             "variable " + variableName +
                                                 " doesn't exist, can't associate attribute " +
                                                 name + ", in call to DefineAttribute");
    }

    const std::string globalName = helper::GlobalName(name, variableName, separator);

    auto itExistingAttribute = m_Attributes.find(globalName);
    if (itExistingAttribute != m_Attributes.end())
    {
        if (itExistingAttribute->second->m_Type == helper::GetDataType<T>())
        {
            if (!itExistingAttribute->second->Equals(static_cast<const void *>(array), elements))
            {
                Attribute<T> &a = static_cast<Attribute<T> &>(*itExistingAttribute->second);
                a.Modify(array, elements);
                void *Data = &a.m_DataSingleValue;
                if (a.m_DataArray.size() != 0)
                    Data = a.m_DataArray.data();
                for (auto &e : m_Engines)
                {
                    e.second->NotifyEngineAttribute(globalName, &a, Data);
                }
            }
        }
        else
        {
            helper::Throw<std::invalid_argument>(
                "Core", "IO", "DefineAttribute",
                "modifiable attribute " + globalName + " has been defined with type " +
                    ToString(itExistingAttribute->second->m_Type) + ". Type cannot be changed to " +
                    ToString(helper::GetDataType<T>()));
        }
        return static_cast<Attribute<T> &>(*itExistingAttribute->second);
    }
    else
    {
        auto itAttributePair =
            m_Attributes.emplace(globalName, std::unique_ptr<AttributeBase>(new Attribute<T>(
                                                 globalName, array, elements, allowModification)));
        Attribute<T> &a = static_cast<Attribute<T> &>(*itAttributePair.first->second);
        void *Data = (void *)array;
        for (auto &e : m_Engines)
        {
            e.second->NotifyEngineAttribute(globalName, &a, Data);
        }
        return static_cast<Attribute<T> &>(*itAttributePair.first->second);
    }
}

template <class T>
Attribute<T> *IO::InquireAttribute(const std::string &name, const std::string &variableName,
                                   const std::string separator) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::InquireAttribute");
    const std::string globalName = helper::GlobalName(name, variableName, separator);
    auto itAttribute = m_Attributes.find(globalName);

    if (itAttribute == m_Attributes.end())
    {
        return nullptr;
    }

    if (itAttribute->second->m_Type != helper::GetDataType<T>())
    {
        return nullptr;
    }

    return static_cast<Attribute<T> *>(itAttribute->second.get());
}

// PRIVATE

template <class T>
Params IO::GetVariableInfo(const std::string &variableName, const std::set<std::string> &keys)
{
    Params info;
    // keys input are case insensitive
    const std::set<std::string> keysLC = helper::LowerCase(keys);

    // return empty map if only "name" key is requested
    if (keys.size() == 1 && keysLC.count("name") == 1)
    {
        return info;
    }

    Variable<T> &variable = *InquireVariable<T>(variableName);

    if (keys.empty() || keysLC.count("type") == 1)
    {
        info["Type"] = ToString(variable.m_Type);
    }

    if (keys.empty() || keysLC.count("availablestepscount") == 1)
    {
        info["AvailableStepsCount"] = helper::ValueToString(variable.m_AvailableStepsCount);
    }

    if (keys.empty() || keysLC.count("shape") == 1)
    {
        // expensive function
        info["Shape"] = helper::VectorToCSV(variable.Shape());
    }

    if (keys.empty() || keysLC.count("singlevalue") == 1)
    {
        const std::string isSingleValue = variable.m_SingleValue ? "true" : "false";
        info["SingleValue"] = isSingleValue;
    }

    if (keys.empty() || (keysLC.count("min") == 1 && keysLC.count("max") == 1))
    {
        if (TypeHasMinMax(helper::GetDataType<T>()))
        {
            const auto pairMinMax = variable.MinMax();
            info["Min"] = helper::ValueToString(pairMinMax.first);
            info["Max"] = helper::ValueToString(pairMinMax.second);
        }
    }
    else if (keysLC.count("min") == 1)
    {
        info["Min"] = helper::ValueToString(variable.Min());
    }
    else if (keysLC.count("max") == 1)
    {
        info["Max"] = helper::ValueToString(variable.Min());
    }
    return info;
}

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_IO_TCC_ */
