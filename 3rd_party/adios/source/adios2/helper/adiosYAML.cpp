/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosYAML.h basic YAML parsing functionality for ADIOS config file schema
 *
 *  Created on: Oct 24, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adiosYAML.h"

#include "adios2/helper/adiosString.h"

#include <yaml-cpp/yaml.h>

namespace adios2
{
namespace helper
{

namespace
{

YAML::Node YAMLNode(const std::string nodeName, const YAML::Node &upperNode,
                    const std::string &hint, const bool isMandatory,
                    const YAML::NodeType::value nodeType)
{
    const YAML::Node node = upperNode[nodeName];

    if (isMandatory && !node)
    {
        helper::Throw<std::invalid_argument>(
            "Helper", "adiosYAML", "YAMLNode",
            "no " + nodeName + " node found, (is your node key lower case?), " + hint);
    }
    if (node && node.Type() != nodeType)
    {
        helper::Throw<std::invalid_argument>("Helper", "adiosYAML", "YAMLNode",
                                             "node " + nodeName +
                                                 " is the wrong type, review adios2 "
                                                 "config YAML specs for the node, " +
                                                 hint);
    }
    return node;
}

Params YAMLNodeMapToParams(const YAML::Node &node, const std::string &hint)
{
    Params parameters;
    for (auto itParam = node.begin(); itParam != node.end(); ++itParam)
    {
        const std::string key = itParam->first.as<std::string>();
        const std::string value = itParam->second.as<std::string>();
        auto it = parameters.emplace(key, value);
        if (!it.second)
        {
            helper::Throw<std::invalid_argument>(
                "Helper", "adiosYAML", "YAMLNodeMapToParams",
                "found duplicated key : " + key + ", keys must be unique in a YAML node, " + hint);
        }
    }
    return parameters;
}

constexpr bool isMandatory = true;
constexpr bool isNotMandatory = false;

inline void FixHomePath(std::string &path, std::string &homePath)
{
    if (!path.empty() && path[0] == '~')
    {
        path = homePath + path.substr(1);
    }
}

/*std::string NodeType(const YAML::Node &node)
{
    switch (node.Type())
    {
    case YAML::NodeType::Null:
        return "Null";
    case YAML::NodeType::Scalar:
        return "Scalar";
    case YAML::NodeType::Sequence:
        return "Sequence";
    case YAML::NodeType::Map:
        return "Map";
    case YAML::NodeType::Undefined:
        return "Undefined";
    }
    return "NoIdeaWhatThisIs";
}*/

template <class T>
void SetOption(T &value, const std::string nodeName, const YAML::Node &upperNode,
               const std::string &hint)
{
    auto node = YAMLNode(nodeName, upperNode, hint, isNotMandatory, YAML::NodeType::Scalar);
    if (node)
    {
        value = node.as<T>();
    }
}

} // end empty  namespace

void IOVariableYAML(const YAML::Node &variableMap, core::IO &currentIO, const std::string &hint)
{
    const YAML::Node &variableNameScalar =
        YAMLNode("Variable", variableMap, hint, isMandatory, YAML::NodeType::Scalar);
    const std::string variableName = variableNameScalar.as<std::string>();

    const YAML::Node operationsSequence =
        YAMLNode("Operations", variableMap, hint, isNotMandatory, YAML::NodeType::Sequence);

    if (operationsSequence)
    {
        // loop through each transport node
        const std::string errorMessage =
            " in operations node from variable " + variableName + ", " + hint;

        for (auto it = operationsSequence.begin(); it != operationsSequence.end(); ++it)
        {
            const YAML::Node typeScalar =
                YAMLNode("Type", *it, errorMessage, isMandatory, YAML::NodeType::Scalar);

            Params parameters = YAMLNodeMapToParams(*it, hint);
            const std::string operatorType = EraseKey<std::string>("Type", parameters);

            currentIO.m_VarOpsPlaceholder[variableName].emplace_back(operatorType, parameters);
        }
    }
}

void IOYAML(core::ADIOS &adios, const YAML::Node &ioMap, core::IO &io, const std::string &hint)
{
    // Engine parameters
    const YAML::Node engineMap = YAMLNode("Engine", ioMap, hint, false, YAML::NodeType::Map);

    if (engineMap)
    {
        Params parameters = YAMLNodeMapToParams(engineMap, hint);
        auto itType = parameters.find("Type");
        if (itType != parameters.end())
        {
            const std::string type = EraseKey<std::string>("Type", parameters);
            io.SetEngine(type);
        }
        io.SetParameters(parameters);
    }

    // Variables
    const YAML::Node variablesSequence =
        YAMLNode("Variables", ioMap, hint, false, YAML::NodeType::Sequence);

    if (variablesSequence)
    {
        // loop through each variable node
        for (const YAML::Node &variableMap : variablesSequence)
        {
            IOVariableYAML(variableMap, io, hint);
        }
    }

    // Transports
    const YAML::Node transportsSequence =
        YAMLNode("Transports", ioMap, hint, false, YAML::NodeType::Sequence);

    if (transportsSequence)
    {
        // loop through each transport node
        for (auto it = transportsSequence.begin(); it != transportsSequence.end(); ++it)
        {
            YAMLNode("Type", *it, " in transport node " + hint, isMandatory,
                     YAML::NodeType::Scalar);

            Params parameters = YAMLNodeMapToParams(*it, hint);
            const std::string type = EraseKey<std::string>("Type", parameters);

            io.AddTransport(type, parameters);
        }
    }
}

void ParseConfigYAMLIO(core::ADIOS &adios, const std::string &configFileYAML,
                       const std::string &configFileContents, core::IO &io)
{
    const std::string hint =
        "when parsing config file " + configFileYAML + " in call to ADIOS constructor";

    // the following copy is needed because YAML::Load modifies configFileContents
    const std::string configFileContentsCopy = configFileContents;
    const YAML::Node document = YAML::Load(configFileContentsCopy);

    if (!document)
    {
        helper::Throw<std::invalid_argument>(
            "Helper", "adiosYAML", "ParseConfigYAML",
            "parser error in file " + configFileYAML +
                " invalid format check with any YAML editor if format is "
                "ill-formed, " +
                hint);
    }

    for (auto itNode = document.begin(); itNode != document.end(); ++itNode)
    {
        const YAML::Node ioScalar =
            YAMLNode("IO", *itNode, hint, isNotMandatory, YAML::NodeType::Scalar);
        if (ioScalar)
        {
            const std::string ioName = ioScalar.as<std::string>();
            if (ioName == io.m_Name)
            {
                IOYAML(adios, *itNode, io, hint);
                return;
            }
        }
    }
}

std::string ParseConfigYAML(core::ADIOS &adios, const std::string &configFileYAML,
                            std::map<std::string, core::IO> &ios)
{
    const std::string hint =
        "when parsing config file " + configFileYAML + " in call to ADIOS constructor";

    const std::string configFileContents = adios.GetComm().BroadcastFile(configFileYAML, hint);
    // the following copy is needed because YAML::Load modifies configFileContents
    const std::string configFileContentsCopy = configFileContents;
    const YAML::Node document = YAML::Load(configFileContentsCopy);

    if (!document)
    {
        helper::Throw<std::invalid_argument>(
            "Helper", "adiosYAML", "ParseConfigYAML",
            "parser error in file " + configFileYAML +
                " invalid format check with any YAML editor if format is "
                "ill-formed, " +
                hint);
    }

    for (auto itNode = document.begin(); itNode != document.end(); ++itNode)
    {
        const YAML::Node ioScalar =
            YAMLNode("IO", *itNode, hint, isNotMandatory, YAML::NodeType::Scalar);
        if (ioScalar)
        {
            const std::string ioName = ioScalar.as<std::string>();
            // Build the IO object
            auto itCurrentIO =
                ios.emplace(std::piecewise_construct, std::forward_as_tuple(ioName),
                            std::forward_as_tuple(adios, ioName, true, adios.m_HostLanguage));
            core::IO &currentIO = itCurrentIO.first->second;
            IOYAML(adios, *itNode, currentIO, hint);
        }
    }
    return configFileContents;
}

void ParseUserOptionsFile(Comm &comm, const std::string &configFileYAML, UserOptions &options,
                          std::string &homePath)
{
    const std::string hint =
        "when parsing user config file " + configFileYAML + " in call to ADIOS constructor";

    const std::string configFileContents = comm.BroadcastFile(configFileYAML, hint);

    const YAML::Node document = YAML::Load(configFileContents);
    if (!document)
    {
        helper::Throw<std::invalid_argument>(
            "Helper", "adiosUserOptions", "ParseUserOptionsFile",
            "parser error in file " + configFileYAML +
                " invalid format. Check with any YAML editor if format is ill-formed, " + hint);
    }

    /*
     * This code section below determines what options we recognize at all from the
     * ~/.config/adios2/adios2.yaml file
     */
    {
        UserOptions::General &opts = options.general;
        const YAML::Node general =
            YAMLNode("General", document, hint, isNotMandatory, YAML::NodeType::Map);
        if (general)
        {
            SetOption(opts.verbose, "verbose", general, hint);
        }
    }

    {
        UserOptions::Campaign &opts = options.campaign;
        const YAML::Node campaign =
            YAMLNode("Campaign", document, hint, isNotMandatory, YAML::NodeType::Map);
        if (campaign)
        {
            SetOption(opts.verbose, "verbose", campaign, hint);
            SetOption(opts.active, "active", campaign, hint);
            SetOption(opts.hostname, "hostname", campaign, hint);
            SetOption(opts.campaignstorepath, "campaignstorepath", campaign, hint);
            FixHomePath(opts.campaignstorepath, homePath);
            SetOption(opts.cachepath, "cachepath", campaign, hint);
            FixHomePath(opts.cachepath, homePath);
        }
    }

    {
        UserOptions::SST &opts = options.sst;
        const YAML::Node sst = YAMLNode("SST", document, hint, isNotMandatory, YAML::NodeType::Map);
        if (sst)
        {
            SetOption(opts.verbose, "verbose", sst, hint);
        }
    }
}

} // end namespace helper
} // end namespace adios2
