/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * AttributeBase.h : base class for Attribute<T> class, allows RTTI at read time
 *
 *  Created on: Aug 1, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_ATTRIBUTEBASE_H_
#define ADIOS2_CORE_ATTRIBUTEBASE_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <string>
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"

namespace adios2
{
namespace core
{

class AttributeBase
{

public:
    const std::string m_Name;
    const DataType m_Type;
    size_t m_Elements;
    bool m_IsSingleValue;
    const bool m_AllowModification;

    /**
     * Single value constructor used by Attribute<T> derived class
     * @param name
     * @param type
     */
    AttributeBase(const std::string &name, const DataType type, const bool allowModification);

    /**
     * Array constructor used by Attribute<T> derived class
     * @param name
     * @param type
     * @param elements
     */
    AttributeBase(const std::string &name, const DataType type, const size_t elements,
                  const bool allowModification);

    virtual ~AttributeBase() = default;

    Params GetInfo() const noexcept;

    /**
     * Compare Attribute's current value with input
     * @return true if they equal, false otherwise
     */
    bool Equals(const void *values, const size_t elements) const noexcept;

private:
    virtual std::string DoGetInfoValue() const noexcept = 0;
    virtual bool DoEqual(const void *values, const size_t elements) const noexcept = 0;
};

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_ATTRIBUTEBASE_H_ */
