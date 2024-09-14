/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_engine.h : Engine public API  C bindings
 *
 *  Created on: Nov 8, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_C_C_ADIOS2_C_ENGINE_H_
#define ADIOS2_BINDINGS_C_C_ADIOS2_C_ENGINE_H_

#include "adios2_c_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Return engine name string and length without '\0\ character
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param name output, string without trailing '\0', NULL or preallocated
 * buffer
 * @param size output, engine_type size without '\0'
 * @param engine handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_engine_name(char *name, size_t *size, const adios2_engine *engine);

/**
 * Return engine type string and length without '\0\ character
 * For safe use, call this function first with NULL name parameter
 * to get the size, then preallocate the buffer (with room for '\0'
 * if desired), then call the function again with the buffer.
 * Then '\0' terminate it if desired.
 * @param type output, string without trailing '\0', NULL or preallocated
 * buffer
 * @param size output, engine_type size without '\0'
 * @param engine handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_engine_get_type(char *type, size_t *size, const adios2_engine *engine);

/**
 * @brief Return the engine's Open mode
 * @param mode output, adios2_mode parameter used in adios2_open()
 * @param engine handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_engine_openmode(adios2_mode *mode, const adios2_engine *engine);

/**
 * @brief Begin a logical adios2 step stream
 * Check each engine documentation for MPI collective/non-collective
 * behavior.
 * @param engine handler
 * @param mode see enum adios2_step_mode in adios2_c_types.h for options,
 * read is the common use case
 * @param timeout_seconds provide a time out in Engine opened in read mode
 * @param status output from enum adios2_step_status in adios2_c_types.h
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_begin_step(adios2_engine *engine, const adios2_step_mode mode,
                               const float timeout_seconds, adios2_step_status *status);
/**
 * Inspect current logical step
 * @param current_step output
 * @param engine input handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_current_step(size_t *current_step, const adios2_engine *engine);

/**
 * Inspect current between step status
 * @param between_step_pairs output boolean
 * @param engine input handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_between_step_pairs(size_t *between_step_pairs, const adios2_engine *engine);

/**
 * Inspect total number of available steps, use for file engines in read mode
 * only
 * @param steps output available steps in engine
 * @param engine input handler
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_steps(size_t *steps, const adios2_engine *engine);

//***************** PUT *****************
/**
 * Put data associated with a Variable in an engine, used for engines with
 * adios2_mode_write at adios2_open
 * @param engine handler for a particular engine where data will be put
 * @param variable contains variable metadata information
 * @param data user data to be associated with a variable, must be the same type
 * passed to adios2_define_variable
 * @param launch mode launch policy
 * <pre>
 * 		adios2_mode_deferred: lazy evaulation, do not use data until
 * first adios2_perform_puts, adios2_end_step, or adios2_close. This is the
 * preferred way.
 * 		adios_mode_sync, data is consumed by the engine and can be
 * reused immediately. Special case, only use if necessary.
 * </pre>
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_put(adios2_engine *engine, adios2_variable *variable, const void *data,
                        const adios2_mode launch);

/**
 * Put data associated with a Variable in an engine, used for engines with
 * adios2_mode_write at adios2_open. This is the name string version
 * @param engine handler for a particular engine where data will be put
 * @param variable_name variable with this name must exists in adios2_io that
 * opened the engine handler (1st parameter)
 * @param data user data to be associated with a variable, must be the same type
 * passed to adios2_define_variable
 * @param launch mode launch policy
 * <pre>
 * 		adios2_mode_deferred: lazy evaulation, do not use data until
 * first adios2_perform_puts, adios2_end_step, or adios2_close. This is the
 * preferred way.
 * 		adios_mode_sync, data is consumed by the engine and can be
 * reused immediately. Special case, only use if necessary.
 * </pre>
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_put_by_name(adios2_engine *engine, const char *variable_name, const void *data,
                                const adios2_mode launch);

/**
 * Performs all the adios2_put and adios2_put_by_name called with mode
 * adios2_mode_deferred, up to this point, by copying user data into internal
 * ADIOS buffers. User data can be reused after this point.
 * @param engine handler for a particular engine where data will be put
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_perform_puts(adios2_engine *engine);

/**
 * Write array data to disk.  This may relieve memory pressure by clearing ADIOS
 * buffers.  It is a collective call. User data can be reused after this point.
 * @param engine handler for a particular engine where data will be put
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_perform_data_write(adios2_engine *engine);

//***************** GET *****************
/**
 * Gets data associated with a Variable from an engine, used for engines with
 * adios2_mode_read at adios2_open.
 * This is the name string version
 * @param engine handler for a particular engine where data will be put
 * @param variable handler must exists in adios2_io that
 * opened the engine handler (1st parameter). Typically from
 * adios2_inquire_variable
 * @param data user data to be associated with a variable, must be the same type
 * passed to adios2_define_variable. Must be pre-allocated for the required
 * variable selection.
 * @param launch mode launch policy
 * <pre>
 * 		adios2_mode_deferred: lazy evaluation, do not use data until
 * first adios2_perform_puts, adios2_end_step, or adios2_close. This is the
 * preferred way.
 * 		adios_mode_sync: data is populated by the engine and can be
 * reused immediately. Special case, only use if necessary.
 * </pre>
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_get(adios2_engine *engine, adios2_variable *variable, void *data,
                        const adios2_mode launch);

/**
 * Gets data associated with a Variable from an engine, used for engines with
 * adios2_mode_read at adios2_open.
 * This is the name string version
 * @param engine handler for a particular engine where data will be put
 * @param variable_name variable with this name must exists in adios2_io that
 * opened the engine handler (1st parameter).
 * @param data user data to be associated with a variable, must be the same type
 * passed to adios2_define_variable. Must be pre-allocated for the required
 * variable selection.
 * @param launch mode launch policy
 * <pre>
 * 		adios2_mode_deferred: lazy evaluation, do not use data until
 * first adios2_perform_puts, adios2_end_step, or adios2_close. This is the
 * preferred way.
 * 		adios_mode_sync, data is populated by the engine and can be
 * reused
 * immediately. Special case, only use if necessary.
 * </pre>
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_get_by_name(adios2_engine *engine, const char *variable_name, void *data,
                                const adios2_mode launch);

/**
 * Performs all the adios2_get and adios2_get_by_name called with mode
 * adios2_mode_deferred up to this point by getting the data from the Engine.
 * User data can be reused after this point.
 * @param engine handler for a particular engine where data will be obtained
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_perform_gets(adios2_engine *engine);

/**
 * Terminates interaction with current step. By default puts/gets data to/from
 * all transports
 * Check each engine documentation for MPI collective/non-collective behavior.
 * @param engine handler executing IO tasks
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_end_step(adios2_engine *engine);

/**
 * Explicit engine buffer flush to transports
 * @param engine input
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_flush(adios2_engine *engine);

/**
 * Explicit engine buffer flush to transport index
 * @param engine input
 * @param transport_index index to be flushed
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_flush_by_index(adios2_engine *engine, const int transport_index);

/**
 * Close all transports in adios2_Engine. Call is required to close system
 * resources.
 * MPI Collective, calls MPI_Comm_free for duplicated communicator at Open
 * @param engine handler containing all transports to
 * be closed. NOTE: engines NEVER become NULL after this function is called.
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_close(adios2_engine *engine);

/**
 * Close a particular transport from the index returned by adios2_add_transport
 * @param engine handler containing all transports to
 * be closed. NOTE: engines NEVER become NULL due to this function.
 * @param transport_index handler from adios2_add_transport
 * @return adios2_error 0: success, see enum adios2_error for errors
 */
adios2_error adios2_close_by_index(adios2_engine *engine, const int transport_index);

/**
 * Promise that no more definitions or changes to defined variables will
 * occur. Useful information if called before the first adios2_end_step() of
 * an output Engine, as it will know that the definitions are complete and
 * constant for the entire lifetime of the output and may optimize metadata
 * handling.
 * @param engine handler
 */
adios2_error adios2_lock_writer_definitions(adios2_engine *engine);

/**
 * Promise that the reader data selections of are fixed and
 * will not change in future timesteps. This information, provided
 * before the EndStep() representing a fixed read pattern, may be
 * utilized by the input Engine to optimize data flow.
 * @param engine handler
 */
adios2_error adios2_lock_reader_selections(adios2_engine *engine);

/**
 * Get the list of blocks for a variable in a given step.
 * In Streaming mode, step is unused, always the current step is processed.
 * @return Newly allocated adios2_varinfo structure, NULL pointer if step does
 * not exist. The memory must be freed by the adios2_free_blockinfo function
 */
adios2_varinfo *adios2_inquire_blockinfo(adios2_engine *engine, adios2_variable *variable,
                                         const size_t step);
/**
 * free adios2_varinfo structure
 * @param data_blocks
 * @return void
 */
void adios2_free_blockinfo(adios2_varinfo *data_blocks);

#ifdef __cplusplus
} // end extern C
#endif

#endif /* ADIOS2_BINDINGS_C_C_ADIOS2_C_ENGINE_H_ */
