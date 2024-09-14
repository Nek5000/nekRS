/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_attribute.h :
 *
 *  Created on: Jun 11, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_C_C_ADIOS2_C_ATTRIBUTE_H_
#define ADIOS2_BINDINGS_C_C_ADIOS2_C_ATTRIBUTE_H_

#include "adios2_c_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Retrieve attribute name
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param name output, string without trailing '\0', NULL or preallocated buffer
 * @param size output, name size without '\0'
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_name(char *name, size_t *size, const adios2_attribute *attribute);

/**
 * Retrieve attribute type
 * @param type
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_type(adios2_type *type, const adios2_attribute *attribute);

/**
 * Retrieve attribute type in string form "char", "unsigned long", etc.
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param type output, string without trailing '\0', NULL or preallocated buffer
 * @param size output, type size without '\0'
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_type_string(char *type, size_t *size,
                                          const adios2_attribute *attribute);

/**
 * Checks if attribute is a single value or an array
 * @param result output, adios2_true: single value, adios2_false: array
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_is_value(adios2_bool *result, const adios2_attribute *attribute);

/**
 * Returns the number of elements (as in C++ STL size() function) if attribute
 * is a 1D array. If single value returns 1
 * @param size output, number of elements in attribute
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_size(size_t *size, const adios2_attribute *attribute);

/**
 * Retrieve attribute data pointer (read-only)
 * @param data output attribute values, must be pre-allocated
 * @param size data size
 * @param attribute handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_attribute_data(void *data, size_t *size, const adios2_attribute *attribute);

#ifdef __cplusplus
} // end extern C
#endif

#endif /* ADIOS2_BINDINGS_C_C_ADIOS2_C_ATTRIBUTE_H_ */
