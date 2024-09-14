/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_io.h
 *
 *  Created on: Nov 8, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_C_C_ADIOS2_C_IO_H_
#define ADIOS2_BINDINGS_C_C_ADIOS2_C_IO_H_

#include "adios2_c_types.h"

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Check if io exists in a config file passed to the adios handler that
 * created this io
 * @param result output adios2_true=1: in config file, adios2_false=0: not in
 * config file
 * @param io handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_in_config_file(adios2_bool *result, const adios2_io *io);

/**
 * @brief Set the engine type for current io handler
 * @param io handler
 * @param engine_type predefined engine type, default is bpfile
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_engine(adios2_io *io, const char *engine_type);

/**
 * @brief Set several parameters at once.
 * @param io handler
 * @param string parameters in the format  "param1=value1 , param2 = value2"
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_parameters(adios2_io *io, const char *parameters);

/**
 * @brief Set a single parameter. Overwrites value if key exists
 * @param io handler
 * @param key parameter key
 * @param value parameter value
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_parameter(adios2_io *io, const char *key, const char *value);

/**
 * Return IO parameter value string and length without '\0\ character
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param value output
 * @param size output, value size without '\0'
 * @param io input handler
 * @param key input parameter key, if not found return size = 0 and value is
 * untouched
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_get_parameter(char *value, size_t *size, const adios2_io *io, const char *key);

/**
 * @brief Clear all parameters.
 * @param io handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_clear_parameters(adios2_io *io);

/**
 * @brief Add a transport to current io handler. Must be supported by current
 * engine type.
 * @param transport_index handler used for setting transport parameters or at
 * adios2_close
 * @param io handler
 * @param type must be a supported transport type for a particular Engine.
 *             CAN'T use the keywords "Transport" or "transport"
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_add_transport(size_t *transport_index, adios2_io *io, const char *type);

/**
 * @brief Set a single parameter to an existing transport identified
 * with a transport_index handler from adios2_add_transport.
 * Overwrites existing parameter with the same key.
 * @param io handler
 * @param transport_index handler from adios2_add_transport
 * @param key parameter key
 * @param value parameter value
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_set_transport_parameter(adios2_io *io, const size_t transport_index,
                                            const char *key, const char *value);

/**
 * @brief Define a variable within io
 * @param io handler that owns the variable
 * @param name unique variable identifier
 * @param type primitive type from enum adios2_type in adios2_c_types.h
 * @param ndims number of dimensions
 * @param shape global dimension
 * @param start local offset
 * @param count local dimension
 * @param constant_dims adios2_constant_dims_true:: shape, start, count
 * won't change; adios2_constant_dims_false: shape, start, count will change
 * after definition
 * @return success: handler, failure: NULL
 */
adios2_variable *adios2_define_variable(adios2_io *io, const char *name, const adios2_type type,
                                        const size_t ndims, const size_t *shape,
                                        const size_t *start, const size_t *count,
                                        const adios2_constant_dims constant_dims);

#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
/**
 * @brief Define a derived variable within io
 * @param io handler that owns the variable
 * @param name unique variable identifier
 * @param type primitive type from enum adios2_type in adios2_c_types.h
 * @param ndims number of dimensions
 * @param shape global dimension
 * @param start local offset
 * @param count local dimension
 * @param constant_dims adios2_constant_dims_true:: shape, start, count
 * won't change; adios2_constant_dims_false: shape, start, count will change
 * after definition
 * @return success: handler, failure: NULL
 */
adios2_derived_variable *adios2_define_derived_variable(adios2_io *io, const char *name,
                                                        const char *expression,
                                                        const adios2_derived_var_type type);
#endif

/**
 * @brief Retrieve a variable handler within current io handler
 * @param io handler to variable io owner
 * @param name unique variable identifier within io handler
 * @return found: handler, not found: NULL
 */
adios2_variable *adios2_inquire_variable(adios2_io *io, const char *name);

/**
 * Returns an array of variable handlers for all variable present in the io
 * group
 * @param variables output array of variable pointers (pointer to an
 * adios2_variable**)
 * @param size output number of variables
 * @param io handler to variables io owner
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_inquire_all_variables(adios2_variable ***variables, size_t *size,
                                          adios2_io *io);

/*
 * list all variables under full_group_name
 */
adios2_error adios2_inquire_group_variables(adios2_variable ***variables,
                                            const char *full_group_name, size_t *size,
                                            adios2_io *io);

/**
 * @brief Define an attribute value inside io
 * @param io handler that owns the attribute
 * @param name unique attribute name inside IO handler
 * @param type primitive type from enum adios2_type in adios2_c_types.h
 * @param value attribute single value
 * @return success: handler, failure: NULL
 */
adios2_attribute *adios2_define_attribute(adios2_io *io, const char *name, const adios2_type type,
                                          const void *value);

/**
 * @brief Define an attribute array inside io
 * @param io handler that owns the attribute
 * @param name unique attribute name inside IO handler
 * @param type primitive type from enum adios2_type in adios2_c_types.h
 * @param data attribute data array
 * @param size number of elements of data array
 * @return success: handler, failure: NULL
 */
adios2_attribute *adios2_define_attribute_array(adios2_io *io, const char *name,
                                                const adios2_type type, const void *data,
                                                const size_t size);

/**
 * Define an attribute single value associated to an existing variable by its
 * name
 * @param io handler that owns the variable and attribute
 * @param name unique attribute name inside a variable in io handler
 * @param type primitive type from enum adios2_type in adios2_c_types.h
 * @param value attribute single value
 * @param variable_name unique variable identifier in io handler. If variable
 * doesn't exist adios2_error is adios2_error_invalid_argument.
 * @param separator hierarchy separator (e.g. "/" in variable_name/name )
 * @return success: handler, failure: NULL
 */
adios2_attribute *adios2_define_variable_attribute(adios2_io *io, const char *name,
                                                   const adios2_type type, const void *value,
                                                   const char *variable_name,
                                                   const char *separator);

/**
 * Define an attribute array associated to an existing variable by its name
 * @param io handler that owns the variable and attribute
 * @param name unique attribute name inside a variable in io handler
 * @param type primitive type from enum adios2_type in adios2_c_types.h
 * @param data attribute data single value or array
 * @param size number of elements of data array
 * @param variable_name unique variable identifier in io handler. If variable
 * doesn't exist adios2_error is true.
 * @param separator hierarchy separator (e.g. "/" in variable/attribute )
 * @return success: handler, failure: NULL
 */
adios2_attribute *adios2_define_variable_attribute_array(adios2_io *io, const char *name,
                                                         const adios2_type type, const void *data,
                                                         const size_t size,
                                                         const char *variable_name,
                                                         const char *separator);

/**
 * Returns a handler to a previously defined attribute by name
 * @param io handler to attribute io owner
 * @param name unique attribute identifier within io handler
 * @return found: handler, not found: NULL
 */
adios2_attribute *adios2_inquire_attribute(adios2_io *io, const char *name);

/**
 * Retrieve a handler to a previously defined attribute associated to a variable
 * @param io handler to attribute and variable io owner
 * @param name unique attribute name inside a variable in io handler
 * @param variable_name name of the variable associate with this attribute
 * @param separator hierarchy separator (e.g. "/" in variable/attribute )
 * @return found: handler, not found: NULL
 */
adios2_attribute *adios2_inquire_variable_attribute(adios2_io *io, const char *name,
                                                    const char *variable_name,
                                                    const char *separator);
/**
 * Returns an array of attribute handlers for all attribute present in the io
 * group
 * @param attributes output array of attribute pointers (pointer to an
 * adios2_attribute**)
 * @param size output number of attributes
 * @param io handler to attributes io owner
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_inquire_all_attributes(adios2_attribute ***attributes, size_t *size,
                                           adios2_io *io);
adios2_error adios2_inquire_group_attributes(adios2_attribute ***attributes,
                                             const char *full_prefix, size_t *size, adios2_io *io);

/**
 * Return a list of list sub group names
 *
 */
adios2_error adios2_inquire_subgroups(char ***subGroupNames, const char *full_prefix, size_t *size,
                                      adios2_io *io);

/**
 * @brief DANGEROUS! Removes a variable identified by name. Might create
 * dangling pointers
 * @param result output adios2_true(1): found and removed variable,
 *                      adios2_false(0): not found, nothing to remove
 * @param io handler variable io owner
 * @param name unique variable name within io handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_remove_variable(adios2_bool *result, adios2_io *io, const char *name);

/**
 * @brief DANGEROUS! Removes all existing variables in current IO object.
 * Might create dangling pointers
 * @param io handler variables io owner
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_remove_all_variables(adios2_io *io);
/**
 * @brief returns an array or c strings for names of available variables
 * Might create dangling pointers
 * @param io handler variables io owner
 * @param length of array of strings
 * @return names of variables as an array of strings
 */
char **adios2_available_variables(adios2_io *io, size_t *size);

/**
 * @brief returns an array or c strings for names of available attributes
 * Might create dangling pointers
 * @param io handler variables io owner
 * @param length of array of strings
 * @return names of variables as an array of strings
 */
char **adios2_available_attributes(adios2_io *io, size_t *size);

/**
 * @brief DANGEROUS! Removes an attribute identified by name. Might create
 * dangling pointers
 * @param result output adios2_true(1): found and removed attribute,
 *                      adios2_false(0): not found, nothing to remove
 * @param io handler attribute io owner
 * @param name unique attribute name within io handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_remove_attribute(adios2_bool *result, adios2_io *io, const char *name);

/**
 * @brief DANGEROUS! Removes all existing attributes in current IO object.
 * Might create dangling pointers
 * @param io handler attributes io owner
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_remove_all_attributes(adios2_io *io);

/**
 * Open an Engine to start heavy-weight input/output operations.
 * In MPI version reuses the communicator from adios2_init or adios2_init_config
 * MPI Collective function as it calls MPI_Comm_dup
 * @param io engine owner
 * @param name unique engine identifier
 * @param mode adios2_mode_write, adios2_mode_read, adios2_mode_append
 * and adios2_mode_readRandomAccess
 * @return success: handler, failure: NULL
 */
adios2_engine *adios2_open(adios2_io *io, const char *name, const adios2_mode mode);

#if ADIOS2_USE_MPI
/**
 * Open an Engine to start heavy-weight input/output operations.
 * MPI Collective function as it calls MPI_Comm_dup
 * @param io engine owner
 * @param name unique engine identifier
 * @param mode adios2_mode_write, adios2_mode_read, adios2_mode_append and
 * adios2_mode_readRandomAccess
 * @param comm communicator other than adios' handler comm. MPI only.
 * @return success: handler, failure: NULL
 */
adios2_engine *adios2_open_new_comm(adios2_io *io, const char *name, const adios2_mode mode,
                                    MPI_Comm comm);
#endif

/**
 * Flushes all engines created with current io handler using adios2_open
 * @param io handler whose engine will be flushed
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_flush_all_engines(adios2_io *io);

/**
 * return engine type string and length without null character
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param engine_type output, string without trailing '\0', NULL or preallocated
 * buffer
 * @param size output, engine_type size without '\0'
 * @param io handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_engine_type(char *engine_type, size_t *size, const adios2_io *io);

adios2_engine *adios2_get_engine(adios2_io *io, const char *name);

#ifdef __cplusplus
} // end extern C
#endif

#endif /* ADIOS2_BINDINGS_C_C_ADIOS2_C_IO_H_ */
