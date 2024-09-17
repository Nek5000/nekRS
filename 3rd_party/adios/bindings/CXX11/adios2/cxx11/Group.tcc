/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Group.tcc :
 *
 *  Created on: August 25, 2020
 *      Author: Dmitry Ganyushin ganyushindi@ornl.gov
 */
#ifndef ADIOS2_BINDINGS_CXX11_CXX11_GROUP_TCC_
#define ADIOS2_BINDINGS_CXX11_CXX11_GROUP_TCC_

#include "Group.h"

#include "adios2/core/Group.h"

namespace adios2
{
template <class T>
Variable<T> Group::InquireVariable(const std::string &name)
{
    helper::CheckForNullptr(m_Group,
                            "for variable name " + name + ", in call to Group::InquireVariable");
    return Variable<T>(m_Group->InquireVariable<typename TypeInfo<T>::IOType>(name));
}

template <class T>
Attribute<T> Group::InquireAttribute(const std::string &name, const std::string &variableName,
                                     const std::string separator)
{
    using IOType = typename TypeInfo<T>::IOType;
    helper::CheckForNullptr(m_Group,
                            "for attribute name " + name + ", in call to IO::InquireAttribute");
    return Attribute<T>(m_Group->InquireAttribute<IOType>(name, variableName, separator));
}

} // end namespace adios2

#endif /* ADIOS2_BINDINGS_CXX11_CXX11_GROUP_TCC_ */
