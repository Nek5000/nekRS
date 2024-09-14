/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Attribute.h : template class that defines typed attributes
 *
 *  Created on: Aug 1, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_ATTRIBUTE_H_
#define ADIOS2_CORE_ATTRIBUTE_H_

#include "adios2/core/AttributeBase.h"

namespace adios2
{
namespace core
{

/** @brief Attributes provide complementary information to IO Variables*/
template <class T>
class Attribute : public AttributeBase
{

public:
    std::vector<T> m_DataArray; ///< holds data for array attributes
    T m_DataSingleValue;        ///< holds data for single value attributes

    /**
     * Copy constructor (enforces zero-padding)
     * @param other
     */
    Attribute(const Attribute<T> &other);

    /**
     * Data array constructor
     * @param name
     * @param data
     * @param elements
     * @param allowModifications
     */
    Attribute(const std::string &name, const T *data, const size_t elements,
              const bool allowModification);

    /**
     * Single value constructor
     * @param name
     * @param data
     * @param elements
     * @param allowModifications
     */
    Attribute(const std::string &name, const T &data, const bool allowModification);

    ~Attribute() = default;

    /**
     * Modification of an existing attribute (array)
     */
    void Modify(const T *data, const size_t elements);

    /**
     * Modification of an existing attribute (single value)
     */
    void Modify(const T &data);

private:
    std::string DoGetInfoValue() const noexcept override;
    bool DoEqual(const void *values, const size_t elements) const noexcept override;
};

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_ATTRIBUTE_H_ */
