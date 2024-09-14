/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_variable.h : exposes some members of the Variable handler
 *
 *  Created on: Nov 10, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_C_C_ADIOS2_C_VARIABLE_H_
#define ADIOS2_BINDINGS_C_C_ADIOS2_C_VARIABLE_H_

#include <stddef.h> //size_t

#include "adios2_c_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Set new shape, care must be taken when reading back the variable for
 * different steps. Only applies to Global arrays.
 * @param variable handler for which new selection will be applied to
 * @param ndims number of dimensions for start and count
 * @param shape new shape dimensions array
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_shape(adios2_variable *variable, const size_t ndims, const size_t *shape);

/**
 * Sets the memory space for all following Puts/Gets
 * to either host (default) or device
 * @param mem memory space where Put/Get buffers are allocated
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_memory_space(adios2_variable *variable, const adios2_memory_space mem);

/**
 * Get the memory space that was set by the application for a given variable
 * @param memory space output, the variable memory space
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_get_memory_space(adios2_memory_space *mem, adios2_variable *variable);

/**
 * Read mode only. Required for reading local variables. For Global Arrays it
 * will Set
 * the appropriate Start and Count Selection for the global array
 * coordinates.
 * @param variable handler for which new selection will be applied to
 * @param block_id variable block index defined at write time. Blocks can be
 * inspected with bpls -D variableName
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_block_selection(adios2_variable *variable, const size_t block_id);

/**
 * Set new start and count dimensions
 * @param variable handler for which new selection will be applied to
 * @param ndims number of dimensions for start and count
 * @param start new start dimensions array
 * @param count new count dimensions array
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_selection(adios2_variable *variable, const size_t ndims,
                                  const size_t *start, const size_t *count);

/**
 * Set the local start (offset) point to the memory pointer passed at Put
 * and the memory local dimensions (count). Used for non-contiguous memory
 * writes and reads (e.g. multidimensional ghost-cells).
 * Currently not working for calls to Get.
 * @param variable handler for which new memory selection will be applied to
 * @param ndims number of dimensions for memory_start and memory_count
 * @param memory_start relative local offset of variable.start to the
 * contiguous memory pointer passed at Put from which data starts. e.g. if
 * variable start = {rank*Ny,0} and there is 1 ghost cell per dimension,
 * then memory_start = {1,1}
 * @param memory_count local dimensions for the contiguous memory pointer
 * passed at adios2_put, e.g. if there is 1 ghost cell per dimension and
 * variable count = {Ny,Nx}, then memory_count = {Ny+2,Nx+2}
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_memory_selection(adios2_variable *variable, const size_t ndims,
                                         const size_t *memory_start, const size_t *memory_count);

/**
 * Set new step selection using step_start and step_count. Used mostly for
 * reading from file-based engines (e.g. bpfile, hdf5)
 * @param variable handler for which new selection will be applied to
 * @param step_start starting step for reading
 * @param step_count number of steps to read from step start
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_step_selection(adios2_variable *variable, const size_t step_start,
                                       const size_t step_count);

/**
 * Retrieve variable name
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param name output, string without trailing '\0', NULL or preallocated buffer
 * @param size output, name size without '\0'
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_name(char *name, size_t *size, const adios2_variable *variable);

/**
 * Retrieve variable type
 * @param type output, from enum adios2_type
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_type(adios2_type *type, const adios2_variable *variable);

/**
 * Retrieve variable type in string form "char", "unsigned long", etc.
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param type output, string without trailing '\0', NULL or preallocated buffer
 * @param size output, type size without '\0'
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_type_string(char *type, size_t *size, const adios2_variable *variable);

/**
 * Retrieve variable shapeid
 * @param shapeid output, from enum adios2_shapeid
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_shapeid(adios2_shapeid *shapeid, const adios2_variable *variable);

/**
 * Retrieve current variable number of dimensions
 * @param ndims output
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_ndims(size_t *ndims, const adios2_variable *variable);

/**
 * Retrieve current variable shape
 * @param shape output, must be pre-allocated with ndims
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_shape(size_t *shape, const adios2_variable *variable);

/**
 * Retrieve current variable shape for a given memory space
 * @param shape output, must be pre-allocated with ndims
 * @param variable handler
 * @param memory space
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_shape_with_memory_space(size_t *shape, const adios2_variable *variable,
                                                     const adios2_memory_space mem);

/**
 * Retrieve current variable start
 * @param start output, single value
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_start(size_t *start, const adios2_variable *variable);

/**
 * Retrieve current variable start
 * @param count output, must be pre-allocated with ndims
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_count(size_t *count, const adios2_variable *variable);

/**
 * Read API, get available steps start from available steps count
 * (e.g. in a file for a variable).
 * @param steps_start output absolute first available step, don't use with
 * adios2_set_step_selection as inputs are relative, use 0 instead.
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_steps_start(size_t *steps_start, const adios2_variable *variable);

/**
 * Read API, get available steps count from available steps count
 * (e.g. in a file for a variable). Not necessarily contiguous.
 * @param steps output available steps, single value
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_steps(size_t *steps, const adios2_variable *variable);

/**
 * Returns the minimum required allocation (in number of elements of a certain
 * type, not bytes) for the current selection
 * @param size number of elements of current type to be allocated by a
 * pointer/vector to read current selection
 * @param variable handler for which data size will be inspected from
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_selection_size(size_t *size, const adios2_variable *variable);

/**
 * Adds an operation to a variable (e.g. compression)
 * @param operation_index output handler to be used with
 adios2_add_operation_param
 * @param variable handler on which operation is applied to
 * @param op handler to adios2_operator associated to current operation
 * @param key parameter key supported by the operation, empty if none
 * @param value parameter value supported by the operation, empty if none
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_add_operation(size_t *operation_index, adios2_variable *variable,
                                  adios2_operator *op, const char *key, const char *value);

/**
 * Adds a parameter to an operation created with adios2_add_operation
 * @param variable handler on which operation is applied to
 * @param operation_id handler returned from adios2_add_operation
 * @param key parameter key supported by the operation
 * @param value parameter value supported by the operation
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_operation_parameter(adios2_variable *variable, const size_t operation_id,
                                            const char *key, const char *value);

/**
 * Removes all current Operations associated with AddOperation.
 * Provides the posibility to apply or not operators on a block basis.
 * @param variable handler on which operation is applied to
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_remove_operations(adios2_variable *variable);

/**
 * Read mode only: return the absolute minimum for current variable
 * @param min output: variable minimum, must be of the same type as the variable
 * handler
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_min(void *min, const adios2_variable *variable);

/**
 * Read mode only: return the absolute maximum for current variable
 * @param max output: variable minimum, must be of the same type as the variable
 * handler
 * @param variable handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_variable_max(void *max, const adios2_variable *variable);

#ifdef __cplusplus
} // end extern C
#endif

#endif /* ADIOS2_BINDINGS_C_C_ADIOS2_C_VARIABLE_H_ */
