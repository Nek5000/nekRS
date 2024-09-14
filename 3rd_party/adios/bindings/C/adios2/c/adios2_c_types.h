/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_types.h : basic adios2 types
 *
 *  Created on: Aug 7, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_C_ADIOS2_C_TYPES_H_
#define ADIOS2_BINDINGS_C_ADIOS2_C_TYPES_H_

#include <stddef.h> // size_t
#include <stdint.h> // SIZE_MAX

#include "adios2/common/ADIOSConfig.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct adios2_adios adios2_adios;
typedef struct adios2_io adios2_io;
typedef struct adios2_variable adios2_variable;
typedef struct adios2_derived_variable adios2_derived_variable;
typedef struct adios2_attribute adios2_attribute;
typedef struct adios2_engine adios2_engine;
typedef struct adios2_operator adios2_operator;

/**
 * @brief adios2_error return types for all ADIOS2 C API functions
 * Based on the library C++ standardized exceptions
 * https://en.cppreference.com/w/cpp/error/exception
 * Each error will issue a more detailed description in the standard error
 * output, stderr
 */
typedef enum
{
    /** success */
    adios2_error_none = 0,

    /**
     * user input error
     */
    adios2_error_invalid_argument = 1,

    /** low-level system error, e.g. system IO error */
    adios2_error_system_error = 2,

    /** runtime errors other than system errors, e.g. memory overflow */
    adios2_error_runtime_error = 3,

    /** any other error exception */
    adios2_error_exception = 4

} adios2_error;

typedef enum
{
    adios2_false = 0,
    adios2_true = 1,
} adios2_bool;

typedef enum
{
    adios2_constant_dims_false = 0,
    adios2_constant_dims_true = 1,
} adios2_constant_dims;

typedef enum
{
    adios2_advance_step_false = 0,
    adios2_advance_step_true = 1,
} adios2_advance_step;

typedef enum
{
    adios2_type_unknown = -1,

    adios2_type_string = 0,
    adios2_type_float = 1,
    adios2_type_double = 2,
    adios2_type_float_complex = 3,
    adios2_type_double_complex = 4,

    adios2_type_int8_t = 5,
    adios2_type_int16_t = 6,
    adios2_type_int32_t = 7,
    adios2_type_int64_t = 8,

    adios2_type_uint8_t = 9,
    adios2_type_uint16_t = 10,
    adios2_type_uint32_t = 11,
    adios2_type_uint64_t = 12,
    adios2_type_long_double = 13 // junmin added
} adios2_type;

typedef enum
{
    adios2_mode_undefined = 0,
    adios2_mode_write = 1,
    adios2_mode_read = 2,
    adios2_mode_append = 3,
    adios2_mode_readRandomAccess = 6,

    adios2_mode_deferred = 4,
    adios2_mode_sync = 5
} adios2_mode;

typedef enum
{
    adios2_step_mode_append = 0,
    adios2_step_mode_update = 1,
    adios2_step_mode_read = 2,
} adios2_step_mode;

typedef enum
{
    adios2_step_status_other_error = -1,
    adios2_step_status_ok = 0,
    adios2_step_status_not_ready = 1,
    adios2_step_status_end_of_stream = 2
} adios2_step_status;

typedef enum
{
    adios2_shapeid_unknown = -1,
    adios2_shapeid_global_value = 0,
    adios2_shapeid_global_array = 1,
    adios2_shapeid_joined_array = 2,
    adios2_shapeid_local_value = 3,
    adios2_shapeid_local_array = 4
} adios2_shapeid;

typedef enum
{
    adios2_arrayordering_rowmajor,
    adios2_arrayordering_columnmajor,
    adios2_arrayordering_auto
} adios2_arrayordering;

#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
/** Type of derived variables */
typedef enum
{
    adios2_derived_var_type_metadata_only = 0,
    adios2_derived_var_type_expression_string = 1,
    adios2_derived_var_type_store_data = 2
} adios2_derived_var_type;
#endif

static const size_t adios2_string_array_element_max_size = 4096;

static const size_t adios2_local_value_dim = SIZE_MAX - 2;

union adios2_PrimitiveStdtypeUnion
{
    int8_t int8;
    int16_t int16;
    int32_t int32;
    int64_t int64;
    uint8_t uint8;
    uint16_t uint16;
    uint32_t uint32;
    uint64_t uint64;
    float f;
    double d;
    long double ld;
    char *str;
};

typedef enum
{
    adios2_memory_space_detect = 0,
    adios2_memory_space_host = 1,
#ifdef ADIOS2_HAVE_GPU_SUPPORT
    adios2_memory_space_gpu = 2,
#endif
} adios2_memory_space;

typedef struct
{
    int WriterID;
    size_t BlockID;
    size_t *Start;
    size_t *Count;
    union adios2_PrimitiveStdtypeUnion MinUnion;
    union adios2_PrimitiveStdtypeUnion MaxUnion;
    union adios2_PrimitiveStdtypeUnion Value;
} adios2_blockinfo;

typedef struct
{
    int Dims;
    size_t *Shape;
    int IsValue;
    int IsReverseDims;
    size_t nblocks;
    adios2_blockinfo *BlocksInfo;
} adios2_varinfo;

#ifdef __cplusplus
} // end extern C
#endif

#endif /* ADIOS2_BINDINGS_C_ADIOS2_C_TYPES_H_ */
