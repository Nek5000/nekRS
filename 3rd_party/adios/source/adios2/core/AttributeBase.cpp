/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * AttributeBase.cpp
 *
 *  Created on: Aug 1, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "AttributeBase.h"

namespace adios2
{
namespace core
{

AttributeBase::AttributeBase(const std::string &name, const DataType type,
                             const bool allowModification)
: m_Name(name), m_Type(type), m_Elements(1), m_IsSingleValue(true),
  m_AllowModification(allowModification)
{
}

AttributeBase::AttributeBase(const std::string &name, const DataType type, const size_t elements,
                             const bool allowModification)
: m_Name(name), m_Type(type), m_Elements(elements), m_IsSingleValue(false),
  m_AllowModification(allowModification)
{
}

Params AttributeBase::GetInfo() const noexcept
{
    Params info;
    info["Type"] = ToString(m_Type);
    info["Elements"] = std::to_string(m_Elements);
    info["Value"] = this->DoGetInfoValue();
    info["Modifiable"] = std::to_string(m_AllowModification);
    return info;
}

bool AttributeBase::Equals(const void *values, const size_t elements) const noexcept
{
    return this->DoEqual(values, elements);
}

} // end namespace core
} // end namespace adios2
