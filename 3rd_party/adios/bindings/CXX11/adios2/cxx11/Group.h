/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Group.h :
 *
 *  Created on: August 25, 2020
 *      Author: Dmitry Ganyushin ganyushindi@ornl.gov
 */
#ifndef ADIOS2_BINDINGS_CXX11_CXX11_GROUP_H_

#define ADIOS2_BINDINGS_CXX11_CXX11_GROUP_H_

#include "Attribute.h"
#include "Variable.h"
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

#include "Group.h"
#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"

namespace adios2
{
class IO;

namespace core
{
class Group; // private implementation
}
class Group
{
    friend class IO;

    Group(core::Group *group);
    core::Group *m_Group = nullptr;

public:
    ~Group();
    /**
     * @brief returns available groups on the path set
     * @param
     * @return vector of strings
     */
    std::vector<std::string> AvailableGroups();
    /**
     * @brief returns available variables on the path set
     * @param
     * @return vector of strings
     */
    std::vector<std::string> AvailableVariables();
    /**
     * @brief returns available attributes on the path set
     * @param
     * @return vector of strings
     */
    std::vector<std::string> AvailableAttributes();
    /**
     * @brief returns the current path
     * @param
     * @return current path as a string
     */
    std::string InquirePath();
    /**
     * @brief set the path, points to a particular node on the tree
     * @param next possible path extension
     */
    void setPath(std::string path);
    /**
     * @brief returns a new group object
     * @param name of the group
     * @return new group object
     */
    Group InquireGroup(std::string group_name);
    /**
     * @brief Gets an existing variable of primitive type by name. A wrapper for
     * the corresponding function of the IO class
     * @param name of variable to be retrieved
     * @return pointer to an existing variable in current IO, nullptr if not
     * found
     */
    template <class T>
    Variable<T> InquireVariable(const std::string &name);
    /**
     * Gets an existing attribute of primitive type by name. A wrapper for
     * the corresponding function of the IO class
     * @param name of attribute to be retrieved
     * @return pointer to an existing attribute in current IO, nullptr if not
     * found
     */
    template <class T>
    Attribute<T> InquireAttribute(const std::string &name, const std::string &variableName = "",
                                  const std::string separator = "/");

    /**
     * Inspects variable type. This function can be used in conjunction with
     * MACROS in an else if (type == adios2::GetType<T>() ) {} loop
     * @param name unique variable name identifier in current IO
     * @return type as in adios2::GetType<T>() (e.g. "double", "float"),
     * empty std::string if variable not found
     */
    DataType VariableType(const std::string &name) const;

    /**
     * Inspects attribute type. This function can be used in conjunction with
     * MACROS in an else if (type == adios2::GetType<T>() ) {} loop
     * @param name unique attribute name identifier in current IO
     * @return type as in adios2::GetType<T>() (e.g. "double", "float"), empty
     * std::string if attribute not found
     */
    DataType AttributeType(const std::string &name) const;
};
}
#endif // ADIOS2_BINDINGS_CXX11_CXX11_GROUP_H_
