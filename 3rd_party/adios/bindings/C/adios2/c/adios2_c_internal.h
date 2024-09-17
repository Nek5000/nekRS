
/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_internal.h: helper functionality for use within the C bindings,
 *                      but not part of the public interface
 *
 *  Created on: Feb 9, 2019
 *      Author: Kai Germaschewski <kai.germaschewski@unh.edu>
 */

#ifndef ADIOS2_BINDINGS_C_C_ADIOS2_C_INTERNAL_H_
#define ADIOS2_BINDINGS_C_C_ADIOS2_C_INTERNAL_H_

#include "adios2_c_types.h"

#include <string>

namespace
{
/**
 * MapAdios2Type
 * maps the C adios2_type enum to the actual C/C++ type,
 * except string, which needs to be handled separately
 */
template <int adios2_type>
struct MapAdios2Type;

adios2_error String2CAPI(const std::string &s, char *buf, size_t *size);

} // anonymous namespace

#include "adios2_c_internal.inl"

#define ADIOS2_FOREACH_C_TYPE_1ARG(MACRO)                                                          \
    MACRO(adios2_type_int8_t)                                                                      \
    MACRO(adios2_type_int16_t)                                                                     \
    MACRO(adios2_type_int32_t)                                                                     \
    MACRO(adios2_type_int64_t)                                                                     \
    MACRO(adios2_type_uint8_t)                                                                     \
    MACRO(adios2_type_uint16_t)                                                                    \
    MACRO(adios2_type_uint32_t)                                                                    \
    MACRO(adios2_type_uint64_t)                                                                    \
    MACRO(adios2_type_float)                                                                       \
    MACRO(adios2_type_double)                                                                      \
    MACRO(adios2_type_long_double)                                                                 \
    MACRO(adios2_type_float_complex)                                                               \
    MACRO(adios2_type_double_complex)

// calls MACRO for all adios2_type attribute types except for adios2_type_string
#define ADIOS2_FOREACH_C_ATTRIBUTE_TYPE_1ARG(MACRO)                                                \
    MACRO(adios2_type_int8_t)                                                                      \
    MACRO(adios2_type_int16_t)                                                                     \
    MACRO(adios2_type_int32_t)                                                                     \
    MACRO(adios2_type_int64_t)                                                                     \
    MACRO(adios2_type_uint8_t)                                                                     \
    MACRO(adios2_type_uint16_t)                                                                    \
    MACRO(adios2_type_uint32_t)                                                                    \
    MACRO(adios2_type_uint64_t)                                                                    \
    MACRO(adios2_type_float)                                                                       \
    MACRO(adios2_type_double)                                                                      \
    MACRO(adios2_type_long_double)                                                                 \
    MACRO(adios2_type_float_complex)                                                               \
    MACRO(adios2_type_double_complex)

#endif /* ADIOS2_BINDINGS_C_C_ADIOS2_C_INTERNAL_H_ */
