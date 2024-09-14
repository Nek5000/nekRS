/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPBase.inl
 *
 *  Created on: Sep 3, 2019
 *      Author: Kai Germaschewski kai.germaschewski@unh.edu
 *              William F Godoy godoywf@ornl.gov
 *
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPBASE_INL_
#define ADIOS2_TOOLKIT_FORMAT_BP_BPBASE_INL_
#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPBASE_H_
#error "Inline file should only be included from its header, never on its own"
#endif

#include "BPBase.h"

namespace adios2
{
namespace format
{

// PROTECTED

#define make_TypeTraits(data_type, TYPE)                                       \
    template <>                                                                \
    struct BPBase::TypeTraits<TYPE>                                            \
    {                                                                          \
        static const DataTypes type_enum = DataTypes::data_type;               \
    };

/* clang-format off */
make_TypeTraits(type_string, std::string)
make_TypeTraits(type_byte, int8_t)
make_TypeTraits(type_short, int16_t)
make_TypeTraits(type_integer, int32_t)
make_TypeTraits(type_long, int64_t)
make_TypeTraits(type_unsigned_byte, uint8_t)
make_TypeTraits(type_unsigned_short, uint16_t)
make_TypeTraits(type_unsigned_integer, uint32_t)
make_TypeTraits(type_unsigned_long, uint64_t)
make_TypeTraits(type_real, float)
make_TypeTraits(type_double, double)
make_TypeTraits(type_long_double, long double)
make_TypeTraits(type_complex, cfloat)
make_TypeTraits(type_double_complex, cdouble)
make_TypeTraits(type_char, char)
/* clang-format on */
#undef make_TypeTraits

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BP_BPBASE_INL_ */
