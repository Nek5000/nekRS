/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO.cpp
 *
 *  Created on: Jan 6, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "IO.h"
#include "IO.tcc"

#include <memory>
#include <mutex>
#include <sstream>
#include <utility> // std::pair

#include "adios2/common/ADIOSMacros.h"

#include "adios2/engine/bp3/BP3Reader.h"
#include "adios2/engine/bp3/BP3Writer.h"
#include "adios2/engine/bp4/BP4Reader.h"
#include "adios2/engine/bp4/BP4Writer.h"
#include "adios2/engine/bp5/BP5Reader.h"
#include "adios2/engine/bp5/BP5Writer.h"
#include "adios2/engine/inline/InlineReader.h"
#include "adios2/engine/inline/InlineWriter.h"
#include "adios2/engine/mhs/MhsReader.h"
#include "adios2/engine/mhs/MhsWriter.h"
#include "adios2/engine/null/NullReader.h"
#include "adios2/engine/null/NullWriter.h"
#include "adios2/engine/plugin/PluginEngine.h"
#include "adios2/engine/skeleton/SkeletonReader.h"
#include "adios2/engine/skeleton/SkeletonWriter.h"

#include "adios2/helper/adiosComm.h"
#include "adios2/helper/adiosFunctions.h" //BuildParametersMap
#include "adios2/helper/adiosString.h"
#include <adios2sys/SystemTools.hxx> // FileIsDirectory()

#ifdef ADIOS2_HAVE_DATAMAN // external dependencies
#include "adios2/engine/dataman/DataManReader.h"
#include "adios2/engine/dataman/DataManWriter.h"
#endif

#ifdef ADIOS2_HAVE_SST // external dependencies
#include "adios2/engine/sst/SstReader.h"
#include "adios2/engine/sst/SstWriter.h"
#endif

#ifdef ADIOS2_HAVE_DAOS // external dependencies
#include "adios2/engine/daos/DaosReader.h"
#include "adios2/engine/daos/DaosWriter.h"
#endif

#ifdef ADIOS2_HAVE_CAMPAIGN // external dependencies
#include "adios2/engine/campaign/CampaignReader.h"
#endif

namespace adios2
{
namespace core
{

IO::EngineFactoryEntry IO_MakeEngine_HDF5();

namespace
{
std::unordered_map<std::string, IO::EngineFactoryEntry> Factory = {
    {"bp3", {IO::MakeEngine<engine::BP3Reader>, IO::MakeEngine<engine::BP3Writer>}},
    {"bp4", {IO::MakeEngine<engine::BP4Reader>, IO::MakeEngine<engine::BP4Writer>}},
    {"bp5", {IO::MakeEngine<engine::BP5Reader>, IO::MakeEngine<engine::BP5Writer>}},
    {"dataman",
#ifdef ADIOS2_HAVE_DATAMAN
     {IO::MakeEngine<engine::DataManReader>, IO::MakeEngine<engine::DataManWriter>}
#else
     IO::NoEngineEntry("ERROR: this version didn't compile with "
                       "DataMan library, can't use DataMan engine\n")
#endif
    },
    {"ssc", IO::NoEngineEntry("ERROR: this version didn't compile with "
                              "SSC library, can't use SSC engine\n")},
    {"mhs",
#ifdef ADIOS2_HAVE_MHS
     {IO::MakeEngine<engine::MhsReader>, IO::MakeEngine<engine::MhsWriter>}
#else
     IO::NoEngineEntry("ERROR: this version didn't compile with "
                       "MHS library, can't use MHS engine\n")
#endif
    },
    {"sst",
#ifdef ADIOS2_HAVE_SST
     {IO::MakeEngine<engine::SstReader>, IO::MakeEngine<engine::SstWriter>}
#else
     IO::NoEngineEntry("ERROR: this version didn't compile with "
                       "Sst library, can't use Sst engine\n")
#endif
    },
    {"daos",
#ifdef ADIOS2_HAVE_DAOS
     {IO::MakeEngine<engine::DaosReader>, IO::MakeEngine<engine::DaosWriter>}
#else
     IO::NoEngineEntry("ERROR: this version didn't compile with "
                       "DAOS library, can't use DAOS engine\n")
#endif
    },
    {"effis",
#ifdef ADIOS2_HAVE_SST
     {IO::MakeEngine<engine::SstReader>, IO::MakeEngine<engine::SstWriter>}
#else
     IO::NoEngineEntry("ERROR: this version didn't compile with "
                       "Sst library, can't use Sst engine\n")
#endif
    },
    {"dataspaces", IO::NoEngineEntry("ERROR: this version didn't compile with "
                                     "DataSpaces library, can't use DataSpaces engine\n")},
    {"hdf5",
#ifdef ADIOS2_HAVE_HDF5
     IO_MakeEngine_HDF5()
#else
     IO::NoEngineEntry("ERROR: this version didn't compile with "
                       "HDF5 library, can't use HDF5 engine\n")
#endif
    },
    {"skeleton", {IO::MakeEngine<engine::SkeletonReader>, IO::MakeEngine<engine::SkeletonWriter>}},
    {"inline", {IO::MakeEngine<engine::InlineReader>, IO::MakeEngine<engine::InlineWriter>}},
    {"null", {IO::MakeEngine<engine::NullReader>, IO::MakeEngine<engine::NullWriter>}},
    {"nullcore",
     {IO::NoEngine("ERROR: nullcore engine does not support read mode"),
      IO::MakeEngine<engine::NullWriter>}},
    {"plugin", {IO::MakeEngine<plugin::PluginEngine>, IO::MakeEngine<plugin::PluginEngine>}},

    {"campaign",
#ifdef ADIOS2_HAVE_CAMPAIGN
     {IO::MakeEngine<engine::CampaignReader>,
      IO::NoEngine("ERROR: campaign engine does not support write mode")}
#else
     IO::NoEngineEntry("ERROR: this version didn't compile with "
                       "support for campaign management, can't use Campaign engine\n")
#endif
    },
};

const std::unordered_map<std::string, bool> ReadRandomAccess_Supported = {
    {"bp3", false},     {"bp4", false},        {"bp5", true},      {"dataman", false},
    {"ssc", false},     {"mhs", false},        {"sst", false},     {"daos", false},
    {"effis", false},   {"dataspaces", false}, {"hdf5", false},    {"skeleton", true},
    {"inline", false},  {"null", true},        {"nullcore", true}, {"plugin", false},
    {"campaign", true},
};

// Synchronize access to the factory in case one thread is
// looking up while another registers additional entries.
std::mutex FactoryMutex;

std::unordered_map<std::string, IO::EngineFactoryEntry>::const_iterator
FactoryLookup(std::string const &name)
{
    std::lock_guard<std::mutex> factoryGuard(FactoryMutex);
    (void)factoryGuard;
    return Factory.find(name);
}

struct ThrowError
{
    std::shared_ptr<Engine> operator()(IO &, const std::string &, const Mode, helper::Comm) const
    {
        helper::Throw<std::invalid_argument>("Core", "IO", "Operator", Err);
        return nullptr;
    }
    std::string Err;
};

} // end anonymous namespace

IO::MakeEngineFunc IO::NoEngine(std::string e) { return ThrowError{e}; }

IO::EngineFactoryEntry IO::NoEngineEntry(std::string e) { return {NoEngine(e), NoEngine(e)}; }

void IO::RegisterEngine(const std::string &engineType, EngineFactoryEntry entry)
{
    std::lock_guard<std::mutex> factoryGuard(FactoryMutex);
    Factory[engineType] = std::move(entry);
}

IO::IO(ADIOS &adios, const std::string name, const bool inConfigFile,
       const std::string hostLanguage)
: m_ADIOS(adios), m_Name(name), m_HostLanguage(hostLanguage), m_InConfigFile(inConfigFile)
{
}

IO::~IO() = default;

void IO::SetEngine(const std::string engineType) noexcept
{
    auto lf_InsertParam = [&](const std::string &key, const std::string &value) {
        m_Parameters.insert(std::pair<std::string, std::string>(key, value));
    };

    /* First step in handling virtual engine names */
    std::string finalEngineType;
    std::string engineTypeLC = engineType;
    std::transform(engineTypeLC.begin(), engineTypeLC.end(), engineTypeLC.begin(), ::tolower);
    if (engineTypeLC == "insituviz" || engineTypeLC == "insituvisualization")
    {
        finalEngineType = "SST";
        lf_InsertParam("FirstTimestepPrecious", "true");
        lf_InsertParam("RendezvousReaderCount", "0");
        lf_InsertParam("QueueLimit", "3");
        lf_InsertParam("QueueFullPolicy", "Discard");
        lf_InsertParam("AlwaysProvideLatestTimestep", "false");
    }
    else if (engineTypeLC == "insituanalysis")
    {
        finalEngineType = "SST";
        lf_InsertParam("FirstTimestepPrecious", "false");
        lf_InsertParam("RendezvousReaderCount", "1");
        lf_InsertParam("QueueLimit", "1");
        lf_InsertParam("QueueFullPolicy", "Block");
        lf_InsertParam("AlwaysProvideLatestTimestep", "false");
    }
    else if (engineTypeLC == "codecoupling")
    {
        finalEngineType = "SST";
        lf_InsertParam("FirstTimestepPrecious", "false");
        lf_InsertParam("RendezvousReaderCount", "1");
        lf_InsertParam("QueueLimit", "1");
        lf_InsertParam("QueueFullPolicy", "Block");
        lf_InsertParam("AlwaysProvideLatestTimestep", "false");
    }
    else if (engineTypeLC == "filestream")
    {
        finalEngineType = "filestream";
        lf_InsertParam("OpenTimeoutSecs", "3600");
        lf_InsertParam("StreamReader", "true");
    }
    /* "file" is handled entirely in IO::Open() as it needs the name */
    else
    {
        finalEngineType = engineType;
    }

    m_EngineType = finalEngineType;
}
void IO::SetIOMode(const IOMode ioMode) { m_IOMode = ioMode; }

void IO::SetParameters(const Params &parameters) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    for (const auto &parameter : parameters)
    {
        m_Parameters[parameter.first] = parameter.second;
    }
}

void IO::SetParameters(const std::string &parameters)
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    adios2::Params parameterMap = adios2::helper::BuildParametersMap(parameters, '=', ',');
    SetParameters(parameterMap);
}

void IO::SetParameter(const std::string key, const std::string value) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    m_Parameters[key] = value;
}

Params &IO::GetParameters() noexcept { return m_Parameters; }

void IO::ClearParameters() noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    m_Parameters.clear();
}

size_t IO::AddTransport(const std::string type, const Params &parameters)
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    Params parametersMap(parameters);

    if (parameters.count("transport") == 1 || parameters.count("Transport") == 1)
    {
        helper::Throw<std::invalid_argument>(
            "Core", "IO", "AddTransport",
            "key Transport (or transport) is not allowed in transport "
            "parameters");
    }

    CheckTransportType(type);

    parametersMap["transport"] = type;
    m_TransportsParameters.push_back(parametersMap);
    return m_TransportsParameters.size() - 1;
}

void IO::SetTransportParameter(const size_t transportIndex, const std::string key,
                               const std::string value)
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    if (transportIndex >= m_TransportsParameters.size())
    {
        helper::Throw<std::invalid_argument>("Core", "IO", "SetTransportParameter",
                                             "transport Index " + std::to_string(transportIndex) +
                                                 " does not exist");
    }

    m_TransportsParameters[transportIndex][key] = value;
}

const VarMap &IO::GetVariables() const noexcept { return m_Variables; }
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
const VarMap &IO::GetDerivedVariables() const noexcept { return m_VariablesDerived; }
#endif

const AttrMap &IO::GetAttributes() const noexcept { return m_Attributes; }

bool IO::InConfigFile() const noexcept { return m_InConfigFile; }

void IO::SetDeclared() noexcept { m_IsDeclared = true; }

void IO::SetArrayOrder(const ArrayOrdering ArrayOrder) noexcept
{
    if (ArrayOrder == ArrayOrdering::Auto)
    {
        if (helper::IsRowMajor(m_HostLanguage))
            m_ArrayOrder = ArrayOrdering::RowMajor;
        else
            m_ArrayOrder = ArrayOrdering::ColumnMajor;
    }
    else
        m_ArrayOrder = ArrayOrder;
}

bool IO::IsDeclared() const noexcept { return m_IsDeclared; }

bool IO::RemoveVariable(const std::string &name) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::RemoveVariable");
    bool isRemoved = false;
    auto itVariable = m_Variables.find(name);
    // variable exists
    if (itVariable != m_Variables.end())
    {
        m_Variables.erase(itVariable);
        isRemoved = true;
    }
    return isRemoved;
}

void IO::RemoveAllVariables() noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::RemoveAllVariables");
    m_Variables.clear();
}

bool IO::RemoveAttribute(const std::string &name) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::RemoveAttribute");
    bool isRemoved = false;
    auto itAttribute = m_Attributes.find(name);
    // attribute exists
    if (itAttribute != m_Attributes.end())
    {
        // first remove the Attribute object
        const DataType type(itAttribute->second->m_Type);

        if (type == DataType::None)
        {
            // nothing to do
        }
        else
        {
            m_Attributes.erase(itAttribute);
            isRemoved = true;
        }
    }

    return isRemoved;
}

void IO::RemoveAllAttributes() noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::RemoveAllAttributes");
    m_Attributes.clear();
}

std::map<std::string, Params> IO::GetAvailableVariables(const std::set<std::string> &keys) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::GetAvailableVariables");

    std::map<std::string, Params> variablesInfo;
    for (const auto &variablePair : m_Variables)
    {
        const std::string variableName = variablePair.first;
        const DataType type = InquireVariableType(variableName);

        if (type == DataType::Struct)
        {
        }
#define declare_template_instantiation(T)                                                          \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        variablesInfo[variableName] = GetVariableInfo<T>(variableName, keys);                      \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
    }

    return variablesInfo;
}

std::map<std::string, Params> IO::GetAvailableAttributes(const std::string &variableName,
                                                         const std::string separator,
                                                         const bool fullNameKeys) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::GetAvailableAttributes");
    std::map<std::string, Params> attributesInfo;

    if (!variableName.empty())
    {
        auto itVariable = m_Variables.find(variableName);
        const DataType type = InquireVariableType(itVariable);

        if (type == DataType::Struct)
        {
        }
        else
        {
            attributesInfo = itVariable->second->GetAttributesInfo(*this, separator, fullNameKeys);
        }
        return attributesInfo;
    }

    // return all attributes if variable name is empty
    for (const auto &attributePair : m_Attributes)
    {
        const std::string &name = attributePair.first;

        if (attributePair.second->m_Type == DataType::Struct)
        {
        }
        else
        {
            attributesInfo[name] = attributePair.second->GetInfo();
        }
    }
    return attributesInfo;
}

DataType IO::InquireVariableType(const std::string &name) const noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    auto itVariable = m_Variables.find(name);
    return InquireVariableType(itVariable);
}

DataType IO::InquireVariableType(const VarMap::const_iterator itVariable) const noexcept
{
    if (itVariable == m_Variables.end())
    {
        return DataType::None;
    }

    const DataType type = itVariable->second->m_Type;

    if (m_ReadStreaming)
    {
        if (type == DataType::Struct)
        {
        }
        else
        {
            if (!itVariable->second->IsValidStep(m_EngineStep + 1))
            {
                return DataType::None;
            }
        }
    }

    return type;
}

DataType IO::InquireAttributeType(const std::string &name, const std::string &variableName,
                                  const std::string separator) const noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    const std::string globalName = helper::GlobalName(name, variableName, separator);

    auto itAttribute = m_Attributes.find(globalName);
    if (itAttribute == m_Attributes.end())
    {
        return DataType::None;
    }

    return itAttribute->second->m_Type;
}

void IO::AddOperation(const std::string &variable, const std::string &operatorType,
                      const Params &parameters) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    m_VarOpsPlaceholder[variable].push_back({operatorType, parameters});
}

Engine &IO::Open(const std::string &name, const Mode mode, helper::Comm comm)
{
    PERFSTUBS_SCOPED_TIMER("IO::Open");
    auto itEngineFound = m_Engines.find(name);
    const bool isEngineFound = (itEngineFound != m_Engines.end());
    bool isEngineActive = false;
    Mode mode_to_use = mode;

    if (isEngineFound)
    {
        if (*itEngineFound->second)
        {
            isEngineActive = true;
        }
    }

    if (isEngineFound)
    {
        if (isEngineActive) // check if active
        {
            helper::Throw<std::invalid_argument>("Core", "IO", "Open",
                                                 "Engine " + name + " is opened twice");
        }
    }

    if (isEngineFound)
    {
        if (!isEngineActive)
        {
            m_Engines.erase(name);
        }
    }

    std::shared_ptr<Engine> engine;
    const bool isDefaultEngine = m_EngineType.empty() ? true : false;
    std::string engineTypeLC = m_EngineType;
    if (!isDefaultEngine)
    {
        std::transform(engineTypeLC.begin(), engineTypeLC.end(), engineTypeLC.begin(), ::tolower);
    }

    /* Second step in handling virtual engines */
    /* BPFile for read needs to use BP5, BP4, or BP3 depending on the file's
     * version
     */
    if ((engineTypeLC == "file" || engineTypeLC == "bpfile" || engineTypeLC == "bp" ||
         isDefaultEngine))
    {
        if (helper::EndsWith(name, ".h5", false))
        {
            engineTypeLC = "hdf5";
        }
        else if (helper::EndsWith(name, ".aca", false))
        {
            engineTypeLC = "campaign";
        }
        else if ((mode_to_use == Mode::Read) || (mode_to_use == Mode::ReadRandomAccess))
        {
            if (adios2sys::SystemTools::FileIsDirectory(name))
            {
                char v = helper::BPVersion(name, comm, m_TransportsParameters);
                engineTypeLC = "bp";
                engineTypeLC.push_back(v);
            }
            else if (adios2sys::SystemTools::FileIsDirectory(name + ".tier0"))
            {
                engineTypeLC = "mhs";
            }
            else
            {
                if (helper::EndsWith(name, ".bp", false))
                {
                    engineTypeLC = "bp3";
                }
                else
                {
                    /* We need to figure out the type of file
                     * from the file itself
                     */
                    if (helper::IsHDF5File(name, *this, comm, m_TransportsParameters))
                    {
                        engineTypeLC = "hdf5";
                    }
                    else
                    {
                        engineTypeLC = "bp3";
                    }
                }
            }
        }
        else
        {
            // File default for writing: BP5
            engineTypeLC = "bp5";
        }
    }

    /* Note: Mismatch between BP4/BP5 writer and FileStream reader is not
       handled if writer has not created the directory yet, when FileStream
       falls back to default (BP4) */
    if (engineTypeLC == "filestream")
    {
        char v = helper::BPVersion(name, comm, m_TransportsParameters);
        engineTypeLC = "bp";
        engineTypeLC.push_back(v);
        // std::cout << "Engine " << engineTypeLC << " selected for FileStream"
        //          << std::endl;
    }

    // For the inline engine, there must be exactly 1 reader, and exactly 1
    // writer.
    if (engineTypeLC == "inline")
    {
        if (mode_to_use == Mode::Append)
        {
            helper::Throw<std::runtime_error>("Core", "IO", "Open",
                                              "Append mode is not supported in the inline engine.");
        }

        // See inline.rst:44
        if (mode_to_use == Mode::Sync)
        {
            helper::Throw<std::runtime_error>("Core", "IO", "Open",
                                              "Sync mode is not supported in the inline engine.");
        }

        if (m_Engines.size() >= 2)
        {
            std::string msg = "Failed to add engine " + name + " to IO \'" + m_Name + "\'. ";
            msg += "An inline engine must have exactly one writer, and one "
                   "reader. ";
            msg += "There are already two engines declared, so no more can be "
                   "added.";
            helper::Throw<std::runtime_error>("Core", "IO", "Open", msg);
        }
        // Now protect against declaration of two writers, or declaration of
        // two readers:
        if (m_Engines.size() == 1)
        {
            auto engine_ptr = m_Engines.begin()->second;
            if (engine_ptr->OpenMode() == mode_to_use)
            {
                std::string msg = "The previously added engine " + engine_ptr->m_Name +
                                  " is already opened in same mode requested for " + name + ". ";
                msg += "The inline engine requires exactly one writer and one "
                       "reader.";
                helper::Throw<std::runtime_error>("Core", "IO", "Open", msg);
            }
        }
    }

    if (mode_to_use == Mode::ReadRandomAccess)
    {
        // older engines don't know about ReadRandomAccess Mode
        auto it = ReadRandomAccess_Supported.find(engineTypeLC);
        if (it != ReadRandomAccess_Supported.end())
        {
            if (!it->second)
            {
                mode_to_use = Mode::Read;
            }
        }
    }

    auto f = FactoryLookup(engineTypeLC);
    if (f != Factory.end())
    {
        if ((mode_to_use == Mode::Read) || (mode_to_use == Mode::ReadRandomAccess))
        {
            engine = f->second.MakeReader(*this, name, mode_to_use, std::move(comm));
        }
        else
        {
            engine = f->second.MakeWriter(*this, name, mode_to_use, std::move(comm));
        }
    }
    else
    {
        helper::Throw<std::invalid_argument>("Core", "IO", "Open",
                                             "Engine type " + m_EngineType + " is not valid");
    }

    auto itEngine = m_Engines.emplace(name, std::move(engine));

    if (!itEngine.second)
    {
        helper::Throw<std::invalid_argument>("Core", "IO", "Open",
                                             "failed to create Engine " + m_EngineType);
    }
    // return a reference
    return *itEngine.first->second.get();
}

Engine &IO::Open(const std::string &name, const Mode mode)
{
    return Open(name, mode, m_ADIOS.GetComm().Duplicate());
}
Group &IO::CreateGroup(char delimiter)
{
    m_Gr = std::make_shared<Group>("", delimiter, *this);
    m_Gr->BuildTree();
    return *m_Gr;
}
Engine &IO::GetEngine(const std::string &name)
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    auto itEngine = m_Engines.find(name);
    if (itEngine == m_Engines.end())
    {
        helper::Throw<std::invalid_argument>("Core", "IO", "GetEngine",
                                             "Engine " + name + " not found");
    }
    // return a reference
    return *itEngine->second.get();
}

void IO::RemoveEngine(const std::string &name)
{
    auto itEngine = m_Engines.find(name);
    if (itEngine != m_Engines.end())
    {
        m_Engines.erase(itEngine);
    }
}

void IO::EnterComputationBlock() noexcept
{
    for (auto &enginePair : m_Engines)
    {
        auto &engine = enginePair.second;
        if (engine->OpenMode() != Mode::Read)
        {
            enginePair.second->EnterComputationBlock();
        }
    }
}

void IO::ExitComputationBlock() noexcept
{
    for (auto &enginePair : m_Engines)
    {
        auto &engine = enginePair.second;
        if (engine->OpenMode() != Mode::Read)
        {
            enginePair.second->ExitComputationBlock();
        }
    }
}

void IO::FlushAll()
{
    PERFSTUBS_SCOPED_TIMER("IO::FlushAll");
    for (auto &enginePair : m_Engines)
    {
        auto &engine = enginePair.second;
        if (engine->OpenMode() != Mode::Read)
        {
            enginePair.second->Flush();
        }
    }
}

void IO::ResetVariablesStepSelection(const bool zeroStart, const std::string hint)
{
    PERFSTUBS_SCOPED_TIMER("IO::other");
    for (auto itVariable = m_Variables.begin(); itVariable != m_Variables.end(); ++itVariable)
    {
        const DataType type = InquireVariableType(itVariable);

        if (type == DataType::None)
        {
            continue;
        }

        if (type == DataType::Struct)
        {
        }
        else
        {
            VariableBase &variable = *itVariable->second;
            variable.CheckRandomAccessConflict(hint);
            variable.ResetStepsSelection(zeroStart);
            variable.m_RandomAccess = false;
        }
    }
}

void IO::SetPrefixedNames(const bool isStep) noexcept
{
    const std::set<std::string> attributes = helper::KeysToSet(m_Attributes);
    const std::set<std::string> variables = helper::KeysToSet(m_Variables);

    for (auto itVariable = m_Variables.begin(); itVariable != m_Variables.end(); ++itVariable)
    {
        // if for each step (BP4), check if variable type is not empty
        // (means variable exist in that step)
        const DataType type = isStep ? InquireVariableType(itVariable) : itVariable->second->m_Type;

        if (type == DataType::None)
        {
            continue;
        }

        if (type == DataType::Struct)
        {
        }
        else
        {
            VariableBase &variable = *itVariable->second;
            variable.m_PrefixedVariables = helper::PrefixMatches(variable.m_Name, variables);
            variable.m_PrefixedAttributes = helper::PrefixMatches(variable.m_Name, attributes);
        }
    }

    m_IsPrefixedNames = true;
}

// PRIVATE
void IO::CheckAttributeCommon(const std::string &name) const
{
    auto itAttribute = m_Attributes.find(name);
    if (itAttribute != m_Attributes.end())
    {
        helper::Throw<std::invalid_argument>("Core", "IO", "CheckAttributeCommon",
                                             "Attribute " + name + " exists in IO " + m_Name +
                                                 ", in call to DefineAttribute");
    }
}

void IO::CheckTransportType(const std::string type) const
{
    if (type.empty() || type.find("=") != type.npos)
    {
        helper::Throw<std::invalid_argument>(
            "Core", "IO", "CheckTransportType",
            "wrong first argument " + type +
                ", must be a single word for a supported transport type, in "
                "call to IO AddTransport");
    }
}

#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
VariableDerived &IO::DefineDerivedVariable(const std::string &name, const std::string &exp_string,
                                           const DerivedVarType varType)
{
    PERFSTUBS_SCOPED_TIMER("IO::DefineDerivedVariable");

    {
        auto itVariable = m_VariablesDerived.find(name);
        if (itVariable != m_VariablesDerived.end())
        {
            helper::Throw<std::invalid_argument>("Core", "IO", "DefineDerivedVariable",
                                                 "derived variable " + name +
                                                     " already defined in IO " + m_Name);
        }
        else
        {
            auto itVariable = m_Variables.find(name);
            if (itVariable != m_Variables.end())
            {
                helper::Throw<std::invalid_argument>(
                    "Core", "IO", "DefineDerivedVariable",
                    "derived variable " + name +
                        " trying to use an already defined variable name in IO " + m_Name);
            }
        }
    }

    derived::Expression derived_exp(exp_string);
    std::vector<std::string> var_list = derived_exp.VariableNameList();
    DataType expressionType = DataType::None;
    bool isConstant = true;
    std::map<std::string, std::tuple<Dims, Dims, Dims>> name_to_dims;
    // check correctness for the variable names and types within the expression
    for (auto var_name : var_list)
    {
        auto itVariable = m_Variables.find(var_name);
        if (itVariable == m_Variables.end())
            helper::Throw<std::invalid_argument>("Core", "IO", "DefineDerivedVariable",
                                                 "using undefine variable " + var_name +
                                                     " in defining the derived variable " + name);
        DataType var_type = InquireVariableType(var_name);
        if (expressionType == DataType::None)
            expressionType = var_type;
        if (expressionType != var_type)
            helper::Throw<std::invalid_argument>("Core", "IO", "DefineDerivedVariable",
                                                 "all variables within a derived variable "
                                                 " must have the same type ");
        if ((itVariable->second)->IsConstantDims() == false)
            isConstant = false;
        name_to_dims.insert({var_name,
                             {(itVariable->second)->m_Start, (itVariable->second)->m_Count,
                              (itVariable->second)->m_Shape}});
    }
    // std::cout << "Derived variable " << name << ": PASS : variables exist and have the same type"
    //          << std::endl;
    // set the initial shape of the expression and check correcness
    derived_exp.SetDims(name_to_dims);
    // std::cout << "Derived variable " << name << ": PASS : initial variable dimensions are valid"
    //          << std::endl;

    // create derived variable with the expression
    auto itVariablePair = m_VariablesDerived.emplace(
        name, std::unique_ptr<VariableBase>(
                  new VariableDerived(name, derived_exp, expressionType, isConstant, varType)));
    VariableDerived &variable = static_cast<VariableDerived &>(*itVariablePair.first->second);

    // check IO placeholder for variable operations
    auto itOperations = m_VarOpsPlaceholder.find(name);
    if (itOperations != m_VarOpsPlaceholder.end())
    {
        // allow to apply an operation only for derived variables that save the data
        if (varType != DerivedVarType::StoreData)
            helper::Throw<std::invalid_argument>(
                "Core", "IO", "DefineDerivedVariable",
                "Operators for derived variables can only be applied "
                " for DerivedVarType::StoreData types.");
        variable.m_Operations.reserve(itOperations->second.size());
        for (auto &operation : itOperations->second)
        {
            variable.AddOperation(operation.first, operation.second);
        }
    }
    return variable;
}
#endif

StructDefinition &IO::DefineStruct(const std::string &name, const size_t size)
{
    return m_ADIOS.m_StructDefinitions.emplace(name, StructDefinition(name, size))->second;
}

VariableStruct &IO::DefineStructVariable(const std::string &name, StructDefinition &def,
                                         const Dims &shape, const Dims &start, const Dims &count,
                                         const bool constantDims)
{

    PERFSTUBS_SCOPED_TIMER("IO::DefineStructVariable");

    {
        auto itVariable = m_Variables.find(name);
        if (itVariable != m_Variables.end())
        {
            helper::Throw<std::invalid_argument>("Core", "IO", "DefineStructVariable",
                                                 "variable " + name + " already defined in IO " +
                                                     m_Name);
        }
    }

    auto itVariablePair =
        m_Variables.emplace(name, std::unique_ptr<VariableBase>(new VariableStruct(
                                      name, def, shape, start, count, constantDims)));

    VariableStruct &variable = static_cast<VariableStruct &>(*itVariablePair.first->second);

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

    def.Freeze();

    return variable;
}

VariableStruct *IO::InquireStructVariable(const std::string &name) noexcept
{
    PERFSTUBS_SCOPED_TIMER("IO::InquireStructVariable");

    if (m_Variables.empty())
    {
        for (auto &e : m_Engines)
        {
            e.second->NotifyEngineNoVarsQuery();
        }
        return nullptr;
    }

    auto itVariable = m_Variables.find(name);

    if (itVariable == m_Variables.end())
    {
        return nullptr;
    }

    if (itVariable->second->m_Type != DataType::Struct)
    {
        return nullptr;
    }

    VariableStruct *variable = static_cast<VariableStruct *>(itVariable->second.get());
    if (m_ReadStreaming)
    {
        if (!variable->IsValidStep(m_EngineStep + 1))
        {
            return nullptr;
        }
    }
    return variable;
}

VariableStruct *IO::InquireStructVariable(const std::string &name, const StructDefinition &def,
                                          const bool allowReorganize) noexcept
{
    auto ret = InquireStructVariable(name);
    if (ret == nullptr)
    {
        return nullptr;
    }

    if (ret->m_WriteStructDefinition->Fields() != def.Fields())
    {
        return nullptr;
    }

    for (size_t i = 0; i < def.Fields(); ++i)
    {
        if (ret->m_WriteStructDefinition->Name(i) != def.Name(i))
        {
            return nullptr;
        }
        if (ret->m_WriteStructDefinition->Offset(i) != def.Offset(i) && !allowReorganize)
        {
            return nullptr;
        }
        if (ret->m_WriteStructDefinition->Type(i) != def.Type(i))
        {
            return nullptr;
        }
        if (ret->m_WriteStructDefinition->ElementCount(i) != def.ElementCount(i))
        {
            return nullptr;
        }
    }

    return ret;
}

// Explicitly instantiate the necessary public template implementations
#define define_template_instantiation(T)                                                           \
    template Variable<T> &IO::DefineVariable<T>(const std::string &, const Dims &, const Dims &,   \
                                                const Dims &, const bool);                         \
    template Variable<T> *IO::InquireVariable<T>(const std::string &) noexcept;

ADIOS2_FOREACH_STDTYPE_1ARG(define_template_instantiation)
#undef define_template_instatiation

#define declare_template_instantiation(T)                                                          \
    template Attribute<T> &IO::DefineAttribute<T>(const std::string &, const T *, const size_t,    \
                                                  const std::string &, const std::string,          \
                                                  const bool);                                     \
    template Attribute<T> &IO::DefineAttribute<T>(                                                 \
        const std::string &, const T &, const std::string &, const std::string, const bool);       \
    template Attribute<T> *IO::InquireAttribute<T>(const std::string &, const std::string &,       \
                                                   const std::string) noexcept;

ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace core
} // end namespace adios2
