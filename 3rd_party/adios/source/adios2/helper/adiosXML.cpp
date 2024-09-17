/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosXML.cpp
 *
 *  Created on: May 17, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 *              Chuck Atkins chuck.atkins@kitware.com
 */

#include "adiosXML.h"

/// \cond EXCLUDE_FROM_DOXYGEN
#include <algorithm> //std::transform
#include <stdexcept> //std::invalid_argument
/// \endcond

#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosXMLUtil.h"

#include <pugixml.hpp>

namespace adios2
{
namespace helper
{

std::string FileContents(core::ADIOS &adios, const std::string &configXML)
{
    const std::string configFileContents(adios.GetComm().BroadcastFile(
        configXML, "when parsing configXML file, in call to the ADIOS constructor"));

    if (configFileContents.empty())
    {
        helper::Throw<std::invalid_argument>("Helper", "AdiosXML", "ParseConfigXML",
                                             "empty config xml file");
    }
    return configFileContents;
}

void OperatorXML(core::ADIOS &adios, const pugi::xml_node &operatorNode, const std::string &hint)
{
    const std::unique_ptr<pugi::xml_attribute> name =
        helper::XMLAttribute("name", operatorNode, hint);

    const std::unique_ptr<pugi::xml_attribute> type =
        helper::XMLAttribute("type", operatorNode, hint);

    std::string typeLowerCase = std::string(type->value());
    std::transform(typeLowerCase.begin(), typeLowerCase.end(), typeLowerCase.begin(), ::tolower);

    const Params parameters = helper::XMLGetParameters(operatorNode, hint);

    adios.DefineOperator(name->value(), typeLowerCase, parameters);
}

void IOVariableXML(const pugi::xml_node &variableNode, core::IO &io, const std::string &hint,
                   std::unordered_map<std::string, std::pair<std::string, Params>> &operators)
{
    const std::string variableName =
        std::string(helper::XMLAttribute("name", variableNode, hint)->value());

    for (const pugi::xml_node &operation : variableNode.children("operation"))
    {
        const std::unique_ptr<pugi::xml_attribute> opName =
            helper::XMLAttribute("operator", operation, hint, false);

        const std::unique_ptr<pugi::xml_attribute> opType =
            helper::XMLAttribute("type", operation, hint, false);

        if (*opName && *opType)
        {
            helper::Throw<std::invalid_argument>(
                "Helper", "AdiosXML", "ParseConfigXML",
                "operator (" + std::string(opName->value()) + ") and type (" +
                    std::string(opType->value()) +
                    ") attributes can't coexist in <operation> element "
                    "inside <variable name=\"" +
                    variableName + "\"> element");
        }

        if (!*opName && !*opType)
        {
            helper::Throw<std::invalid_argument>("Helper", "AdiosXML", "ParseConfigXML",
                                                 "<operation> element "
                                                 "inside <variable name=\"" +
                                                     variableName +
                                                     "\"> element requires either operator "
                                                     "(existing) or type (supported) attribute");
        }

        std::string type;
        Params params;

        if (*opName)
        {
            auto itOperator = operators.find(std::string(opName->value()));
            if (itOperator == operators.end())
            {
                helper::Throw<std::invalid_argument>("Helper", "AdiosXML", "ParseConfigXML",
                                                     "operator " + std::string(opName->value()) +
                                                         " not previously defined, from variable " +
                                                         variableName + " inside io " + io.m_Name);
            }
            type = itOperator->second.first;
            params = itOperator->second.second;
        }

        if (*opType)
        {
            type = std::string(opType->value());
        }

        for (const auto &p : helper::XMLGetParameters(operation, hint))
        {
            params[p.first] = p.second;
        }

        io.m_VarOpsPlaceholder[variableName].emplace_back(type, params);
    }
}

void IOXML(core::ADIOS &adios, const pugi::xml_node &ioNode, core::IO &io, const std::string &hint,
           std::unordered_map<std::string, std::pair<std::string, Params>> &operators)
{
    // must be unique per io
    const std::unique_ptr<pugi::xml_node> engine =
        helper::XMLNode("engine", ioNode, hint, false, true);

    if (*engine)
    {
        const std::unique_ptr<pugi::xml_attribute> type =
            helper::XMLAttribute("type", *engine, hint);
        io.SetEngine(type->value());

        const Params parameters = helper::XMLGetParameters(*engine, hint);
        io.SetParameters(parameters);
    }

    for (const pugi::xml_node &transport : ioNode.children("transport"))
    {
        const std::unique_ptr<pugi::xml_attribute> type =
            helper::XMLAttribute("type", transport, hint);

        const Params parameters = helper::XMLGetParameters(transport, hint);
        io.AddTransport(type->value(), parameters);
    }

    for (const pugi::xml_node &variable : ioNode.children("variable"))
    {
        IOVariableXML(variable, io, hint, operators);
    }
}

void ParseConfigXMLIO(core::ADIOS &adios, const std::string &configFileXML,
                      const std::string &configFileContents, core::IO &io,
                      std::unordered_map<std::string, std::pair<std::string, Params>> &operators)
{
    const std::string hint("for config file " + configFileXML + " in call to ADIOS constructor");

    // the following copy is needed because pugi::xml_document modifies configFileContents
    const std::string configFileContentsCopy = configFileContents;
    const std::unique_ptr<pugi::xml_document> document =
        helper::XMLDocument(configFileContentsCopy, hint);

    // must be unique
    const std::unique_ptr<pugi::xml_node> config =
        helper::XMLNode("adios-config", *document, hint, true);

    for (const pugi::xml_node &ioNode : config->children("io"))
    {
        const std::unique_ptr<pugi::xml_attribute> ioName =
            helper::XMLAttribute("name", ioNode, hint);
        if (io.m_Name == ioName->value())
        {
            IOXML(adios, ioNode, io, hint, operators);
            return;
        }
    }
}

std::string
ParseConfigXML(core::ADIOS &adios, const std::string &configFileXML,
               std::map<std::string, core::IO> &ios,
               std::unordered_map<std::string, std::pair<std::string, Params>> &operators)
{
    const std::string hint("for config file " + configFileXML + " in call to ADIOS constructor");

    const std::string configFileContents = FileContents(adios, configFileXML);
    // the following copy is needed because pugi::xml_document modifies configFileContents
    const std::string configFileContentsCopy = configFileContents;
    const std::unique_ptr<pugi::xml_document> document =
        helper::XMLDocument(configFileContentsCopy, hint);

    // must be unique
    const std::unique_ptr<pugi::xml_node> config =
        helper::XMLNode("adios-config", *document, hint, true);

    for (const pugi::xml_node &opNode : config->children("operator"))
    {
        OperatorXML(adios, opNode, hint);
    }

    for (const pugi::xml_node &ioNode : config->children("io"))
    {
        const std::unique_ptr<pugi::xml_attribute> ioName =
            helper::XMLAttribute("name", ioNode, hint);
        // Build the IO object
        auto itCurrentIO =
            ios.emplace(std::piecewise_construct, std::forward_as_tuple(ioName->value()),
                        std::forward_as_tuple(adios, ioName->value(), true, adios.m_HostLanguage));
        core::IO &currentIO = itCurrentIO.first->second;
        IOXML(adios, ioNode, currentIO, hint, operators);
    }
    return configFileContents;
}

} // end namespace helper
} // end namespace adios2
