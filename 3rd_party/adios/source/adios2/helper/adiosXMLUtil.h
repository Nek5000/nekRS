/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosXMLUnit.h: basic XML parsing functionality using pugixml
 *                 forward declares pugi object to keep the pugi dependency
 *                 private in adios2
 *
 *  Created on: April 28, 2020
 *      Author: William F Godoy godoywf@ornl.gov
 */
#ifndef ADIOS2_HELPER_ADIOSXMLUTIL_H_
#define ADIOS2_HELPER_ADIOSXMLUTIL_H_

#include <memory> //std::unique_ptr
#include <string>

#include "adios2/common/ADIOSTypes.h"

// forward declaring pugi objects for private implementation
namespace pugi
{
class xml_attribute;
class xml_node;
class xml_document;
};

namespace adios2
{
namespace helper
{

/**
 * Create a pugi document object from a string with XML contents
 * @param xmlContents input xml document contents in single string
 * @param hint used for exception information
 * @return pugi::xml_document must be unique_ptr to keep pugixml linked
 * privately
 */
std::unique_ptr<pugi::xml_document> XMLDocument(const std::string &xmlContents,
                                                const std::string hint);

/**
 * Extract a XML node from a document (first level node)
 * @param nodeName input node name in document
 * @param xmlDocument input XML document
 * @param hint used for exception information
 * @param isMandatory true: throws exception if node is required and not found,
 * false: not mandatory
 * @param isUnique true: throws exception if node is not unique, false: not
 * unique (many nodes)
 * @return pugi::xml_node must be unique_ptr to keep pugixml linked privately
 */
std::unique_ptr<pugi::xml_node> XMLNode(const std::string nodeName,
                                        const pugi::xml_document &xmlDocument,
                                        const std::string hint, const bool isMandatory = true,
                                        const bool isUnique = false);

/**
 * Extract a XML node from another node (immediate upper-level node)
 * @param nodeName input node name in document
 * @param xmlDocument input XML document to parse
 * @param hint used for exception information
 * @param isMandatory true: throws exception if node is required and not found,
 * false: not mandatory
 * @param isUnique true: throws exception if node is not unique, false: not
 * unique (many nodes)
 * @return pugi::xml_node must be unique_ptr to keep pugixml linked privately
 */
std::unique_ptr<pugi::xml_node> XMLNode(const std::string nodeName, const pugi::xml_node &upperNode,
                                        const std::string hint, const bool isMandatory = true,
                                        const bool isUnique = false);

/**
 * Extract a XML attribute from a node
 * @param attributeName input attribute name in node
 * @param node input node to parse
 * @param hint used for exception information
 * @param isMandatory true: throws exception if attribute is required and not
 * found, false: not mandatory
 * @return pugi::xml_attribute must be unique_ptr to keep pugixml linked
 * privately
 */
std::unique_ptr<pugi::xml_attribute> XMLAttribute(const std::string attributeName,
                                                  const pugi::xml_node &node,
                                                  const std::string hint,
                                                  const bool isMandatory = true);

/**
 * Get Parameters map of key/value strings from a XML node
 * @param node input XML node
 * @param hint used for exception information
 * @param key identifier in the XML schema for parameter key
 * @param value identifier in the XML schema for parameter value
 * @return parameters map
 */
adios2::Params XMLGetParameters(const pugi::xml_node &node, const std::string hint,
                                const std::string xmlKey = "key",
                                const std::string xmlValue = "value");

} // end namespace helper
} // end namespace adios2

#endif /* ADIOS2_HELPER_ADIOSXMLUTIL_H_ */
