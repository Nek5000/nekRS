/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ExampleWritePlugin.cpp
 *
 *  Created on: Jul 5, 2021
 *      Author: Chuck Atkins <chuck.atkins@kitware.com>
 *              Caitlin Ross <caitlin.ross@kitware.com>
 */

#include "ExampleWritePlugin.h"
#include "ExampleWritePlugin.tcc"

#include "adios2/helper/adiosSystem.h"

namespace adios2
{
namespace plugin
{

ExampleWritePlugin::ExampleWritePlugin(core::IO &io, const std::string &name, const Mode mode,
                                       helper::Comm comm)
: PluginEngineInterface(io, name, mode, comm.Duplicate())
{
    Init();
}

ExampleWritePlugin::~ExampleWritePlugin()
{
    m_DataFile.close();
    m_VarFile.close();
}

void ExampleWritePlugin::Init()
{
    std::string dir = "ExamplePlugin";
    auto paramFileNameIt = m_IO.m_Parameters.find("DirName");
    if (paramFileNameIt != m_IO.m_Parameters.end())
    {
        dir = paramFileNameIt->second;
    }
    helper::CreateDirectory(dir);

    std::string fileName = dir + "/data.txt";
    m_DataFile.open(fileName);
    if (!m_DataFile)
    {
        throw std::ios_base::failure("ExampleWritePlugin: Failed to open file " + fileName);
    }

    std::string varfName = dir + "/vars.txt";
    m_VarFile.open(varfName);
    if (!m_VarFile)
    {
        throw std::ios_base::failure("ExampleWritePlugin: Failed to open file " + varfName);
    }
}

StepStatus ExampleWritePlugin::BeginStep(StepMode mode, const float timeoutSeconds)
{
    WriteVarsFromIO();
    return StepStatus::OK;
}

size_t ExampleWritePlugin::CurrentStep() const { return m_CurrentStep; }

void ExampleWritePlugin::EndStep() { m_CurrentStep++; }

#define declare(T)                                                                                 \
    void ExampleWritePlugin::DoPutSync(core::Variable<T> &variable, const T *values)               \
    {                                                                                              \
        WriteArray(variable, values);                                                              \
    }                                                                                              \
    void ExampleWritePlugin::DoPutDeferred(core::Variable<T> &variable, const T *values)           \
    {                                                                                              \
        WriteArray(variable, values);                                                              \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare)
#undef declare

void ExampleWritePlugin::PerformPuts() { WriteVarsFromIO(); }

void ExampleWritePlugin::DoClose(const int transportIndex) {}

void ExampleWritePlugin::WriteVarsFromIO()
{
    const core::VarMap &variables = m_IO.GetVariables();
    for (const auto &vpair : variables)
    {
        const std::string &varName = vpair.first;
        const DataType varType = vpair.second->m_Type;
#define declare_template_instantiation(T)                                                          \
    if (varType == helper::GetDataType<T>())                                                       \
    {                                                                                              \
        core::Variable<T> *v = m_IO.InquireVariable<T>(varName);                                   \
        if (!v)                                                                                    \
            return;                                                                                \
        WriteVariableInfo(*v);                                                                     \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
    }
}

} // end namespace plugin
} // end namespace adios2

extern "C" {

adios2::plugin::ExampleWritePlugin *EngineCreate(adios2::core::IO &io, const std::string &name,
                                                 const adios2::Mode mode, adios2::helper::Comm comm)
{
    return new adios2::plugin::ExampleWritePlugin(io, name, mode, comm.Duplicate());
}

void EngineDestroy(adios2::plugin::ExampleWritePlugin *obj) { delete obj; }
}
