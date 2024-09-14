/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_f2c_operator.cpp :
 *
 *  Created on: Jul 30, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adios2_f2c_common.h"

#ifdef __cplusplus
extern "C" {
#endif

void FC_GLOBAL(adios2_operator_type_f2c,
               ADIOS2_OPERATOR_TYPE_F2C)(char *type, const adios2_operator **op, int *ierr)
{
    size_t sizeC;
    *ierr = static_cast<int>(adios2_operator_type(type, &sizeC, *op));
}

void FC_GLOBAL(adios2_operator_type_length_f2c,
               ADIOS2_OPERATOR_TYPE_LENGTH_F2C)(int *size, const adios2_operator **op, int *ierr)
{
    *size = -1;
    size_t sizeC;
    *ierr = static_cast<int>(adios2_operator_type(nullptr, &sizeC, *op));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *size = static_cast<int>(sizeC);
    }
}

#ifdef __cplusplus
}
#endif
