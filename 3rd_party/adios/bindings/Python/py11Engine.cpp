/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Engine.cpp
 *
 *  Created on: Mar 15, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "py11Engine.h"

#include "adios2/common/ADIOSMacros.h"
#include "adios2/core/Engine.h"
#include "adios2/helper/adiosFunctions.h"

#include <sstream>

#include "py11types.h"

namespace adios2
{
namespace py11
{

Engine::Engine(core::Engine *engine) : m_Engine(engine) {}

Engine::operator bool() const noexcept
{
    if (m_Engine == nullptr)
    {
        return false;
    }

    return *m_Engine ? true : false;
}

StepStatus Engine::BeginStep(const StepMode mode, const float timeoutSeconds)
{
    helper::CheckForNullptr(m_Engine, "in call to Engine::BeginStep");
    return m_Engine->BeginStep(mode, timeoutSeconds);
}

StepStatus Engine::BeginStep()
{
    helper::CheckForNullptr(m_Engine, "in call to Engine::BeginStep");
    return m_Engine->BeginStep();
}

void Engine::Put(Variable variable, const pybind11::array &array, const Mode launch)
{
    helper::CheckForNullptr(m_Engine, "in call to Engine::Put numpy array");
    helper::CheckForNullptr(variable.m_VariableBase,
                            "for variable, in call to Engine::Put numpy array");

    const adios2::DataType type = helper::GetDataTypeFromString(variable.Type());

    if (type == adios2::DataType::Struct)
    {
        // not supported
    }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        m_Engine->Put(*dynamic_cast<core::Variable<T> *>(variable.m_VariableBase),                 \
                      reinterpret_cast<const T *>(array.data()), launch);                          \
    }
    ADIOS2_FOREACH_NUMPY_TYPE_1ARG(declare_type)
#undef declare_type
    else
    {
        throw std::invalid_argument("ERROR: for variable " + variable.Name() +
                                    " numpy array type " + variable.Type() +
                                    " is not supported (found type " + ToString(type) +
                                    ") or "
                                    "is not memory contiguous "
                                    ", in call to Put\n");
    }
}

void Engine::Put(Variable variable, const std::vector<int64_t> &ints, const Mode launch)
{
    helper::CheckForNullptr(m_Engine, "in call to Engine::Put list of ints");
    helper::CheckForNullptr(variable.m_VariableBase,
                            "for variable, in call to Engine::Put list of ints");

    m_Engine->Put(*dynamic_cast<core::Variable<int64_t> *>(variable.m_VariableBase),
                  reinterpret_cast<const int64_t *>(ints.data()), launch);
}

void Engine::Put(Variable variable, const std::vector<double> &floats, const Mode launch)
{
    helper::CheckForNullptr(m_Engine, "in call to Engine::Put list of floats");
    helper::CheckForNullptr(variable.m_VariableBase,
                            "for variable, in call to Engine::Put list of floats");

    m_Engine->Put(*dynamic_cast<core::Variable<double> *>(variable.m_VariableBase),
                  reinterpret_cast<const double *>(floats.data()), launch);
}

void Engine::Put(Variable variable, const std::vector<std::complex<double>> &complexes,
                 const Mode launch)
{
    helper::CheckForNullptr(m_Engine, "in call to Engine::Put list of complexes");
    helper::CheckForNullptr(variable.m_VariableBase,
                            "for variable, in call to Engine::Put list of complexes");
    m_Engine->Put(*dynamic_cast<core::Variable<std::complex<double>> *>(variable.m_VariableBase),
                  reinterpret_cast<const std::complex<double> *>(complexes.data()), launch);
}

void Engine::Put(Variable variable, const std::string &string)
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::Put string");
    helper::CheckForNullptr(variable.m_VariableBase, "for variable, in call to Engine::Put string");

    if (helper::GetDataTypeFromString(variable.Type()) != helper::GetDataType<std::string>())
    {
        throw std::invalid_argument("ERROR: variable " + variable.Name() +
                                    " is not of string type, in call to Engine::Put");
    }

    m_Engine->Put(*dynamic_cast<core::Variable<std::string> *>(variable.m_VariableBase), string,
                  adios2::Mode::Sync);
}

void Engine::PerformPuts()
{
    helper::CheckForNullptr(m_Engine, "in call to PerformPuts");
    m_Engine->PerformPuts();
}

void Engine::PerformDataWrite()
{
    helper::CheckForNullptr(m_Engine, "in call to PerformDataWrite");
    m_Engine->PerformDataWrite();
}

void Engine::Get(Variable variable, pybind11::array &array, const Mode launch)
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::Get a numpy array");
    helper::CheckForNullptr(variable.m_VariableBase,
                            "for variable, in call to Engine::Get a numpy array");

    const adios2::DataType type = helper::GetDataTypeFromString(variable.Type());

    if (type == adios2::DataType::Struct)
    {
        // not supported
    }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        if (!array.dtype().is(pybind11::dtype::of<T>()))                                           \
        {                                                                                          \
            throw std::invalid_argument("In ADIOS2 Get - Type mismatch between Python buffer and " \
                                        "incoming data.");                                         \
        }                                                                                          \
        m_Engine->Get(*dynamic_cast<core::Variable<T> *>(variable.m_VariableBase),                 \
                      reinterpret_cast<T *>(const_cast<void *>(array.data())), launch);            \
    }
    ADIOS2_FOREACH_NUMPY_TYPE_1ARG(declare_type)
#undef declare_type
    else
    {
        throw std::invalid_argument("ERROR: in variable " + variable.Name() + " of type " +
                                    variable.Type() +
                                    ", numpy array type is 1) not supported, 2) a type mismatch or"
                                    "3) is not memory contiguous "
                                    ", in call to Get\n");
    }
}

std::string Engine::Get(Variable variable, const Mode launch)
{
    std::string string;
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::Get a numpy array");
    helper::CheckForNullptr(variable.m_VariableBase,
                            "for variable, in call to Engine::Get a string");

    const adios2::DataType type = helper::GetDataTypeFromString(variable.Type());

    if (type == helper::GetDataType<std::string>())
    {
        m_Engine->Get(*dynamic_cast<core::Variable<std::string> *>(variable.m_VariableBase), string,
                      launch);
    }
    else
    {
        throw std::invalid_argument("ERROR: variable " + variable.Name() + " of type " +
                                    variable.Type() + " is not string, in call to Engine::Get");
    }
    return string;
}
void Engine::PerformGets()
{
    helper::CheckForNullptr(m_Engine, "in call to Engine::PerformGets");
    m_Engine->PerformGets();
}

void Engine::EndStep()
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::EndStep");
    m_Engine->EndStep();
}

bool Engine::BetweenStepPairs() const
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::EndStep");
    return m_Engine->BetweenStepPairs();
}

void Engine::Flush(const int transportIndex)
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::Flush");
    m_Engine->Flush(transportIndex);
}

void Engine::Close(const int transportIndex)
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::Close");
    m_Engine->Close(transportIndex);

    // erase Engine object from IO
    core::IO &io = m_Engine->GetIO();
    const std::string name = m_Engine->m_Name;
    io.RemoveEngine(name);
    m_Engine = nullptr;
}

size_t Engine::CurrentStep() const
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::CurrentStep");
    return m_Engine->CurrentStep();
}

std::string Engine::Name() const
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::Name");
    return m_Engine->m_Name;
}

std::string Engine::Type() const
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::Type");
    return m_Engine->m_EngineType;
}

size_t Engine::Steps() const
{
    helper::CheckForNullptr(m_Engine, "for engine, in call to Engine::Steps");
    return m_Engine->Steps();
}

void Engine::LockWriterDefinitions() const
{
    helper::CheckForNullptr(m_Engine, "in call to Engine::LockWriterDefinitions");
    m_Engine->LockWriterDefinitions();
}

void Engine::LockReaderSelections() const
{
    helper::CheckForNullptr(m_Engine, "in call to Engine::LockReaderSelections");
    m_Engine->LockReaderSelections();
}

std::vector<std::map<std::string, std::string>> Engine::BlocksInfo(std::string &var_name,
                                                                   const size_t step) const
{
    std::vector<std::map<std::string, std::string>> rv;
    auto &varMap = m_Engine->m_IO.GetVariables();
    auto itVariable = varMap.find(var_name);
    if (itVariable == varMap.end())
    {
        return rv;
    }

    // Grab the specified variable object and get its type string
    adios2::DataType var_type = m_Engine->GetIO().InquireVariableType(var_name);

    MinVarInfo *minBlocksInfo = nullptr;

    auto Variable = itVariable->second.get();
    minBlocksInfo = m_Engine->MinBlocksInfo(*Variable, 0);
    if (minBlocksInfo)
    {
        for (auto &info : minBlocksInfo->BlocksInfo)
        {
            std::map<std::string, std::string> info_map;
            std::stringstream start_ss;
            if (info.Start == nullptr)
            {
                start_ss << "0";
            }
            else
            {
                for (size_t i = 0; i < (size_t)minBlocksInfo->Dims; ++i)
                {
                    if (i)
                    {
                        start_ss << ",";
                    }
                    start_ss << (minBlocksInfo->WasLocalValue ? reinterpret_cast<size_t>(info.Start)
                                                              : info.Start[i]);
                }
            }
            info_map["Start"] = start_ss.str();
            std::stringstream count_ss;
            if (info.Count == nullptr)
            {
                count_ss << "0";
            }
            else
            {
                for (size_t i = 0; i < (size_t)minBlocksInfo->Dims; ++i)
                {
                    if (i)
                    {
                        count_ss << ",";
                    }
                    count_ss << (minBlocksInfo->WasLocalValue ? reinterpret_cast<size_t>(info.Count)
                                                              : info.Count[i]);
                }
            }
            info_map["Count"] = count_ss.str();
            info_map["WriterID"] = std::to_string(info.WriterID);
            info_map["BlockID"] = std::to_string(info.BlockID);
            info_map["IsValue"] = minBlocksInfo->IsValue ? "True" : "False";
            std::ostringstream osMax, osMin;
            switch (var_type)
            {
            case DataType::Int8:
                osMax << info.MinMax.MaxUnion.field_int8;
                osMin << info.MinMax.MinUnion.field_int8;
                break;
            case DataType::Int16:
                osMax << info.MinMax.MaxUnion.field_int16;
                osMin << info.MinMax.MinUnion.field_int16;
                break;
            case DataType::Int32:
                osMax << info.MinMax.MaxUnion.field_int32;
                osMin << info.MinMax.MinUnion.field_int32;
                break;
            case DataType::Int64:
                osMax << info.MinMax.MaxUnion.field_int64;
                osMin << info.MinMax.MinUnion.field_int64;
                break;
            case DataType::UInt8:
                osMax << info.MinMax.MaxUnion.field_uint8;
                osMin << info.MinMax.MinUnion.field_uint8;
                break;
            case DataType::UInt16:
                osMax << info.MinMax.MaxUnion.field_uint16;
                osMin << info.MinMax.MinUnion.field_uint16;
                break;
            case DataType::UInt32:
                osMax << info.MinMax.MaxUnion.field_uint32;
                osMin << info.MinMax.MinUnion.field_uint32;
                break;
            case DataType::UInt64:
                osMax << info.MinMax.MaxUnion.field_uint64;
                osMin << info.MinMax.MinUnion.field_uint64;
                break;
            case DataType::Float:
                osMax << info.MinMax.MaxUnion.field_float;
                osMin << info.MinMax.MinUnion.field_float;
                break;
            case DataType::Double:
                osMax << info.MinMax.MaxUnion.field_double;
                osMin << info.MinMax.MinUnion.field_double;
                break;
            case DataType::LongDouble:
                osMax << info.MinMax.MaxUnion.field_ldouble;
                osMin << info.MinMax.MinUnion.field_ldouble;
                break;
            default:
                break;
            }
            info_map["Max"] = osMax.str();
            info_map["Min"] = osMin.str();
            info_map["IsReverseDims"] = minBlocksInfo->IsReverseDims ? "True" : "False";
            rv.push_back(info_map);
        }
        delete minBlocksInfo;
        return rv;
    }
    // Use the macro incantation to call the right instantiation of
    // core::BlocksInfo<>() Note that we are flatting the Dims type items, and
    // returning everything as a dictionary of strings.
    if (false)
    {
    }
#define GET_BLOCKS_INFO(T)                                                                         \
    else if (var_type == helper::GetDataType<T>())                                                 \
    {                                                                                              \
        auto variable = m_Engine->GetIO().InquireVariable<T>(var_name);                            \
        auto infoVec = m_Engine->BlocksInfo<T>(*variable, step);                                   \
        for (auto &info : infoVec)                                                                 \
        {                                                                                          \
            std::map<std::string, std::string> info_map;                                           \
            std::stringstream start_ss;                                                            \
            for (size_t i = 0; i < info.Start.size(); ++i)                                         \
            {                                                                                      \
                if (i != 0)                                                                        \
                    start_ss << ",";                                                               \
                start_ss << info.Start[i];                                                         \
            }                                                                                      \
            info_map["Start"] = start_ss.str();                                                    \
            std::stringstream count_ss;                                                            \
            for (size_t i = 0; i < info.Count.size(); ++i)                                         \
            {                                                                                      \
                if (i != 0)                                                                        \
                    count_ss << ",";                                                               \
                count_ss << info.Count[i];                                                         \
            }                                                                                      \
            info_map["Count"] = count_ss.str();                                                    \
            info_map["WriterID"] = std::to_string(info.WriterID);                                  \
            info_map["BlockID"] = std::to_string(info.BlockID);                                    \
            info_map["IsValue"] = info.IsValue ? "True" : "False";                                 \
            std::ostringstream osMax, osMin;                                                       \
            osMax << info.Max;                                                                     \
            osMin << info.Min;                                                                     \
            info_map["Max"] = osMax.str();                                                         \
            info_map["Min"] = osMin.str();                                                         \
            info_map["IsReverseDims"] = info.IsReverseDims ? "True" : "False";                     \
            rv.push_back(info_map);                                                                \
        }                                                                                          \
    }

    ADIOS2_FOREACH_PYTHON_TYPE_1ARG(GET_BLOCKS_INFO)
#undef GET_BLOCKS_INFO
    else
    {
        throw std::invalid_argument("ERROR: variable " + var_name +
                                    " can't be defined, either type is not "
                                    "supported or is not memory "
                                    "contiguous, in call to DefineVariable\n");
    }

    return rv;
}

} // end namespace py11
} // end namespace adios2
