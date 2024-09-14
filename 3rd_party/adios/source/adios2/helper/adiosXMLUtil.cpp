/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosXML.h basic XML parsing functionality for ADIOS config file schema
 *
 *  Created on: May 17, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adiosXMLUtil.h"
#include "adiosLog.h"

#include <iterator>  // std::distance
#include <stdexcept> //std::invalid_argument

#include <pugixml.hpp>

namespace adios2
{

namespace helper
{

std::unique_ptr<pugi::xml_document> XMLDocument(const std::string &xmlContents,
                                                const std::string hint)
{
    std::unique_ptr<pugi::xml_document> document(new pugi::xml_document);
    auto parse_result =
        document->load_buffer_inplace(const_cast<char *>(xmlContents.data()), xmlContents.size());

    if (!parse_result)
    {
        helper::Throw<std::invalid_argument>(
            "Helper", "adiosXMLUtil", "XMLDocument",
            "parse error in XML string, description: " + std::string(parse_result.description()) +
                ", check with any XML editor if format is ill-formed, " + hint);
    }
    return document;
}

std::unique_ptr<pugi::xml_node> XMLNode(const std::string nodeName,
                                        const pugi::xml_document &xmlDocument,
                                        const std::string hint, const bool isMandatory,
                                        const bool isUnique)
{
    std::unique_ptr<pugi::xml_node> node(new pugi::xml_node(xmlDocument.child(nodeName.c_str())));

    if (isMandatory && !node)
    {
        helper::Throw<std::invalid_argument>("Helper", "adiosXMLUtil", "XMLNode",
                                             "no <" + nodeName + "> element found, " + hint);
    }

    if (isUnique)
    {
        const size_t nodes = std::distance(xmlDocument.children(nodeName.c_str()).begin(),
                                           xmlDocument.children(nodeName.c_str()).end());
        if (nodes > 1)
        {
            helper::Throw<std::invalid_argument>("Helper", "adiosXMLUtil", "XMLNode",
                                                 "XML only one <" + nodeName +
                                                     "> element can exist inside " +
                                                     std::string(xmlDocument.name()) + ", " + hint);
        }
    }
    return node;
}

std::unique_ptr<pugi::xml_node> XMLNode(const std::string nodeName, const pugi::xml_node &upperNode,
                                        const std::string hint, const bool isMandatory,
                                        const bool isUnique)
{
    std::unique_ptr<pugi::xml_node> node(new pugi::xml_node(upperNode.child(nodeName.c_str())));

    if (isMandatory && !node)
    {
        helper::Throw<std::invalid_argument>("Helper", "adiosXMLUtil", "XMLNode",
                                             "no <" + nodeName + "> element found, inside <" +
                                                 std::string(upperNode.name()) + "> element " +
                                                 hint);
    }

    if (isUnique)
    {
        const size_t nodes = std::distance(upperNode.children(nodeName.c_str()).begin(),
                                           upperNode.children(nodeName.c_str()).end());
        if (nodes > 1)
        {
            helper::Throw<std::invalid_argument>(
                "Helper", "adiosXMLUtil", "XMLNode",
                "XML only one <" + nodeName + "> element can exist inside <" +
                    std::string(upperNode.name()) + "> element, " + hint);
        }
    }
    return node;
}

std::unique_ptr<pugi::xml_attribute> XMLAttribute(const std::string attributeName,
                                                  const pugi::xml_node &node,
                                                  const std::string hint, const bool isMandatory)
{
    std::unique_ptr<pugi::xml_attribute> attribute(
        new pugi::xml_attribute(node.attribute(attributeName.c_str())));

    if (isMandatory && !attribute)
    {
        const std::string nodeName(node.name());

        helper::Throw<std::invalid_argument>("Helper", "adiosXMLUtil", "XMLAttribute",
                                             "No attribute " + attributeName + " found on <" +
                                                 nodeName + "> element" + hint);
    }
    return attribute;
}

adios2::Params XMLGetParameters(const pugi::xml_node &node, const std::string hint,
                                const std::string xmlKey, const std::string xmlValue)
{
    const std::string errorMessage("in node " + std::string(node.value()) + ", " + hint);
    Params parameters;

    for (const pugi::xml_node &paramNode : node.children("parameter"))
    {
        const std::unique_ptr<pugi::xml_attribute> key =
            XMLAttribute("key", paramNode, errorMessage);

        const std::unique_ptr<pugi::xml_attribute> value =
            XMLAttribute("value", paramNode, errorMessage);

        parameters.emplace(key->value(), value->value());
    }
    return parameters;
}

} // end namespace helper
} // end namespace adios2
