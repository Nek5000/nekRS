/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Group.cpp :
 *
 *  Created on: August 25, 2020
 *      Author: Dmitry Ganyushin ganyushindi@ornl.gov
 */
#include "Group.h"
#include "Group.tcc"
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "adios2/common/ADIOSMacros.h"
#include "adios2/core/IO.h"
namespace adios2
{
namespace core
{
std::vector<std::string> split(const std::string &s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}
void Group::setPath(std::string path) { currentPath = ADIOS_root + "/" + path; }
void Group::setDelimiter(char delimiter) { groupDelimiter = delimiter; }

Group::Group(std::string path, char delimiter, IO &io)
: currentPath(path), groupDelimiter(delimiter), m_IO(io)
{
    if (mapPtr == nullptr)
    {
        mapPtr = std::shared_ptr<TreeMap>(new TreeMap());
    }
}
// copy constructor
Group::Group(const Group &G)
: currentPath(G.currentPath), groupDelimiter(G.groupDelimiter), m_IO(G.m_IO)
{
    mapPtr = G.mapPtr;
}
Group *Group::InquireGroup(std::string groupName)
{
    if (currentPath.compare("") != 0)
    {
        groupName = currentPath + groupDelimiter + groupName;
    }
    m_Gr = std::make_shared<Group>(groupName, this->groupDelimiter, this->m_IO);
    m_Gr->mapPtr = this->mapPtr;
    return m_Gr.get();
}

void Group::BuildTree()
{
    const core::VarMap &variables = m_IO.GetVariables();
    for (const auto &variablePair : variables)
    {
        std::vector<std::string> tokens = split(variablePair.first, groupDelimiter);
        // Adding artificial root element
        if (tokens[0] == "")
            tokens[0] = ADIOS_root;
        else
            tokens.insert(tokens.begin(), ADIOS_root);
        currentPath = ADIOS_root;
        if (tokens.size() > 1)
        {
            std::string key = tokens[0];
            for (size_t level = 1; level < tokens.size(); level++)
            {
                std::string value = tokens[level];
                // get previous vector
                std::set<std::string> val = mapPtr->treeMap[key];
                // modify it
                val.insert(value);
                mapPtr->treeMap[key] = val;
                key += groupDelimiter + tokens[level];
            }
        }
    }
    const core::AttrMap &attributes = m_IO.GetAttributes();
    for (const auto &attributePair : attributes)
    {
        std::vector<std::string> tokens = split(attributePair.first, groupDelimiter);
        if (tokens.size() > 1)
        {
            std::string key = tokens[0];
            for (size_t level = 1; level < tokens.size(); level++)
            {
                std::string value = tokens[level];
                // get previous vector
                std::set<std::string> val = mapPtr->treeMap[key];
                // modify it
                val.insert(value);
                mapPtr->treeMap[key] = val;
                key += groupDelimiter + tokens[level];
            }
        }
    }
}
std::vector<std::string> Group::AvailableVariables()
{
    // look into map
    std::set<std::string> val = mapPtr->treeMap[currentPath];
    std::vector<std::string> available_variables;
    for (auto const &v : val)
    {
        if (mapPtr->treeMap.find(currentPath + groupDelimiter + v) == mapPtr->treeMap.end())
        {
            const core::VarMap &variables = m_IO.GetVariables();
            std::string variablePath = currentPath + groupDelimiter + v;
            variablePath =
                variablePath.substr(ADIOS_root.size() + 1, variablePath.size() - ADIOS_root.size());
            if (variables.find(variablePath) != variables.end())
            {
                available_variables.push_back(v);
            }
        }
    }

    return available_variables;
}

std::vector<std::string> Group::AvailableAttributes()
{
    // look into map
    std::set<std::string> val = mapPtr->treeMap[currentPath];
    // TODODG check that currentPath exists
    std::vector<std::string> available_attributes;
    for (auto const &v : val)
    {
        if (mapPtr->treeMap.find(currentPath + groupDelimiter + v) == mapPtr->treeMap.end())
        {
            const core::AttrMap &attributes = m_IO.GetAttributes();
            std::string variablePath = currentPath + groupDelimiter + v;
            variablePath =
                variablePath.substr(ADIOS_root.size() + 1, variablePath.size() - ADIOS_root.size());
            if (attributes.find(variablePath) != attributes.end())
            {
                available_attributes.push_back(v);
            }
        }
    }

    return available_attributes;
}

std::vector<std::string> Group::AvailableGroups()
{

    std::vector<std::string> available_groups;
    std::set<std::string> val = mapPtr->treeMap[currentPath];
    for (auto const &v : val)
    {
        if (mapPtr->treeMap.find(currentPath + groupDelimiter + v) != mapPtr->treeMap.end())
            available_groups.push_back(v);
    }
    return available_groups;
}

std::map<std::string, std::set<std::string>> &Group::getTreeMap()
{
    std::map<std::string, std::set<std::string>> &tree = mapPtr->treeMap;
    return tree;
}

std::string Group::InquirePath() { return currentPath; }
Group::~Group() = default;
DataType Group::InquireVariableType(const std::string &name) const noexcept
{

    return m_IO.InquireVariableType(currentPath + groupDelimiter + name);
}

DataType Group::InquireAttributeType(const std::string &name, const std::string &variableName,
                                     const std::string separator) const noexcept
{
    return m_IO.InquireAttributeType(name, variableName, separator);
}
// Explicitly instantiate the necessary public template implementations
#define define_template_instantiation(T)                                                           \
    template Variable<T> *Group::InquireVariable<T>(const std::string &) noexcept;

ADIOS2_FOREACH_STDTYPE_1ARG(define_template_instantiation)
#undef define_template_instatiation

} // end namespace core
} // end namespace adios2
