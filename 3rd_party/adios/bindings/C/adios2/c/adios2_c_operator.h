/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_operator.h :
 *
 *  Created on: Jul 30, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_C_C_ADIOS2_C_OPERATOR_H_
#define ADIOS2_BINDINGS_C_C_ADIOS2_C_OPERATOR_H_

#include "adios2_c_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Retrieve operator type
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param type output, string without trailing '\0', NULL or preallocated buffer
 * @param size output, type size without '\0'
 * @param op operator handler to be inspected
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_operator_type(char *type, size_t *size, const adios2_operator *op);

#ifdef __cplusplus
} // end extern C
#endif

#endif /* ADIOS2_BINDINGS_C_C_ADIOS2_C_OPERATOR_H_ */
