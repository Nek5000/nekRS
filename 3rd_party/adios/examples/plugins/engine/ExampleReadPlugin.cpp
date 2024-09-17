/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ExampleReadPlugin.cpp
 *
 *  Created on: Jul 5, 2021
 *      Author: Caitlin Ross <caitlin.ross@kitware.com>
 */

#include "ExampleReadPlugin.h"
#include "ExampleReadPlugin.tcc"

#include "adios2/helper/adiosType.h"

namespace adios2
{
namespace plugin
{

ExampleReadPlugin::ExampleReadPlugin(core::IO &io, const std::string &name, const Mode mode,
                                     helper::Comm comm)
: PluginEngineInterface(io, name, mode, comm.Duplicate())
{
    Init();
}

ExampleReadPlugin::~ExampleReadPlugin()
{
    m_DataFile.close();
    m_VarFile.close();
}

Dims convertStrToDims(const std::string &str)
{
    Dims dims;
    if (str.size() > 2)
    {
        auto vals = str.substr(1, str.size() - 2);
        std::stringstream ss(vals);
        while (ss.good())
        {
            std::string substr;
            std::getline(ss, substr, ',');
            dims.push_back(std::stoi(substr));
        }
    }
    return dims;
}

void ExampleReadPlugin::Init()
{
    std::string dir = "ExamplePlugin";
    auto paramFileNameIt = m_IO.m_Parameters.find("DirName");
    if (paramFileNameIt != m_IO.m_Parameters.end())
    {
        dir = paramFileNameIt->second;
    }

    std::string fileName = dir + "/data.txt";
    m_DataFile.open(fileName, std::ofstream::in);
    if (!m_DataFile)
    {
        throw std::ios_base::failure("ExampleReadPlugin: Failed to open file " + fileName);
    }

    std::string varfName = dir + "/vars.txt";
    m_VarFile.open(varfName, std::ofstream::in);
    if (!m_VarFile)
    {
        throw std::ios_base::failure("ExampleReadPlugin: Failed to open file " + varfName +
                                     ".vars");
    }

    // get var info
    while (m_VarFile.good())
    {
        std::string name, typeStr, shapeStr, startStr, countStr;
        std::getline(m_VarFile, name, ';');
        std::getline(m_VarFile, typeStr, ';');
        std::getline(m_VarFile, shapeStr, ';');
        std::getline(m_VarFile, startStr, ';');
        std::getline(m_VarFile, countStr);

        auto shape = convertStrToDims(shapeStr);
        auto start = convertStrToDims(startStr);
        auto count = convertStrToDims(countStr);

        const DataType type = helper::GetDataTypeFromString(typeStr);
        if (type == DataType::Struct)
        {
            // not supported
        }
#define declare_template_instantiation(T)                                                          \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        AddVariable<T>(name, shape, start, count);                                                 \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
    }
}

#define declare(T)                                                                                 \
    void ExampleReadPlugin::DoGetSync(core::Variable<T> &variable, T *values)                      \
    {                                                                                              \
        ReadVariable(variable, values);                                                            \
    }                                                                                              \
    void ExampleReadPlugin::DoGetDeferred(core::Variable<T> &variable, T *values)                  \
    {                                                                                              \
        ReadVariable(variable, values);                                                            \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare)
#undef declare

StepStatus ExampleReadPlugin::BeginStep(StepMode mode, const float timeoutSeconds)
{
    return StepStatus::OK;
}

void ExampleReadPlugin::PerformGets() {}

size_t ExampleReadPlugin::CurrentStep() const { return m_CurrentStep; }

void ExampleReadPlugin::EndStep() { m_CurrentStep++; }

void ExampleReadPlugin::DoClose(const int transportIndex) {}

} // end namespace plugin
} // end namespace adios2

extern "C" {

adios2::plugin::ExampleReadPlugin *EngineCreate(adios2::core::IO &io, const std::string &name,
                                                const adios2::Mode mode, adios2::helper::Comm comm)
{
    return new adios2::plugin::ExampleReadPlugin(io, name, mode, comm.Duplicate());
}

void EngineDestroy(adios2::plugin::ExampleReadPlugin *obj) { delete obj; }
}
