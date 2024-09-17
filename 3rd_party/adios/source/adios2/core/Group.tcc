/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Group.tcc template implementations with fix types and specializations
 *
 *  Created on: August 25, 2020
 *      Author: Dmitry Ganyushin ganyushindi@ornl.gov
 */
#ifndef ADIOS2_CORE_GROUP_TCC_
#define ADIOS2_CORE_GROUP_TCC_

#include "Group.h"

namespace adios2
{
namespace core
{

template <class T>
Variable<T> *Group::InquireVariable(const std::string &name) noexcept
{
    std::string variablePath = currentPath + groupDelimiter + name;
    variablePath =
        variablePath.substr(ADIOS_root.size() + 1, variablePath.size() - ADIOS_root.size());
    Variable<T> &variable = *m_IO.InquireVariable<T>(variablePath);
    return &variable;
}

template <class T>
Attribute<T> *Group::InquireAttribute(const std::string &name, const std::string &variableName,
                                      const std::string separator) noexcept
{
    std::string variablePath = currentPath + groupDelimiter + name;
    variablePath =
        variablePath.substr(ADIOS_root.size() + 1, variablePath.size() - ADIOS_root.size());
    Attribute<T> &attribute = m_IO.InquireAttribute<T>(variablePath, variableName, separator);
    return &attribute;
}
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_GROUP_TCC_ */
