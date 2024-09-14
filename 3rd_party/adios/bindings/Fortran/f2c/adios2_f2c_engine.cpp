/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_f2c_engine.cpp
 *
 *  Created on: Nov 8, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adios2_f2c_common.h"

#ifdef __cplusplus
extern "C" {
#endif

void FC_GLOBAL(adios2_begin_step_f2c,
               ADIOS2_BEGIN_STEP_F2C)(adios2_engine **engine, const int *step_mode,
                                      const float *timeout_seconds, int *status, int *ierr)
{
    *status = -1;
    adios2_step_status statusC;

    *ierr = static_cast<int>(adios2_begin_step(*engine, static_cast<adios2_step_mode>(*step_mode),
                                               *timeout_seconds, &statusC));

    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *status = static_cast<int>(statusC);
    }
}

void FC_GLOBAL(adios2_current_step_f2c,
               ADIOS2_CURRENT_STEP_F2C)(int64_t *step, const adios2_engine **engine, int *ierr)
{
    *step = -1;
    size_t stepC;
    *ierr = static_cast<int>(adios2_current_step(&stepC, *engine));

    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *step = static_cast<int64_t>(stepC);
    }
}

void FC_GLOBAL(adios2_steps_f2c, ADIOS2_STEPS_F2C)(int64_t *steps, const adios2_engine **engine,
                                                   int *ierr)
{
    *steps = -1;
    size_t stepsC;
    *ierr = static_cast<int>(adios2_steps(&stepsC, *engine));

    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *steps = static_cast<int64_t>(stepsC);
    }
}

void FC_GLOBAL(adios2_lock_writer_definitions_f2c,
               ADIOS2_LOCK_WRITER_DEFINITIONS_F2C)(adios2_engine **engine, int *ierr)
{
    *ierr = static_cast<int>(adios2_lock_writer_definitions(*engine));
}

void FC_GLOBAL(adios2_lock_reader_selections_f2c,
               ADIOS2_LOCK_READER_SELECTIONS_F2C)(adios2_engine **engine, int *ierr)
{
    *ierr = static_cast<int>(adios2_lock_reader_selections(*engine));
}

// ******** PUTS */
void FC_GLOBAL(adios2_put_f2c, ADIOS2_PUT_F2C)(adios2_engine **engine, adios2_variable **variable,
                                               const void *data, const int *launch, int *ierr)
{
    *ierr =
        static_cast<int>(adios2_put(*engine, *variable, data, static_cast<adios2_mode>(*launch)));
}

void FC_GLOBAL(adios2_put_by_name_f2c, ADIOS2_PUT_BY_NAME_F2C)(adios2_engine **engine,
                                                               const char *name, const void *data,
                                                               const int *launch, int *ierr)
{
    *ierr = static_cast<int>(
        adios2_put_by_name(*engine, name, data, static_cast<adios2_mode>(*launch)));
}

void FC_GLOBAL(adios2_perform_puts_f2c, ADIOS2_PERFORM_PUTS_F2C)(adios2_engine **engine, int *ierr)
{
    *ierr = static_cast<int>(adios2_perform_puts(*engine));
}

void FC_GLOBAL(adios2_perform_data_write_f2c, ADIOS2_PERFORM_DATA_WRITE_F2C)(adios2_engine **engine,
                                                                             int *ierr)
{
    *ierr = static_cast<int>(adios2_perform_data_write(*engine));
}

// ******** GETS */
void FC_GLOBAL(adios2_get_f2c, ADIOS2_get_F2C)(adios2_engine **engine, adios2_variable **variable,
                                               void *data, const int *launch, int *ierr)
{
    *ierr =
        static_cast<int>(adios2_get(*engine, *variable, data, static_cast<adios2_mode>(*launch)));
}

void FC_GLOBAL(adios2_get_by_name_f2c, ADIOS2_get_BY_NAME_F2C)(adios2_engine **engine,
                                                               const char *name, void *data,
                                                               const int *launch, int *ierr)
{
    *ierr = static_cast<int>(
        adios2_get_by_name(*engine, name, data, static_cast<adios2_mode>(*launch)));
}

void FC_GLOBAL(adios2_perform_gets_f2c, ADIOS2_PERFORM_GETS_F2C)(adios2_engine **engine, int *ierr)
{
    *ierr = static_cast<int>(adios2_perform_gets(*engine));
}

void FC_GLOBAL(adios2_end_step_f2c, ADIOS2_END_STEP_F2C)(adios2_engine **engine, int *ierr)
{
    *ierr = static_cast<int>(adios2_end_step(*engine));
}

void FC_GLOBAL(adios2_flush_f2c, ADIOS2_FLUSH_F2C)(adios2_engine **engine, int *ierr)
{
    *ierr = static_cast<int>(adios2_flush(*engine));
}

void FC_GLOBAL(adios2_close_f2c, ADIOS2_CLOSE_F2C)(adios2_engine **engine, int *ierr)
{
    *ierr = static_cast<int>(adios2_close(*engine));
}

void FC_GLOBAL(adios2_engine_get_type_f2c,
               ADIOS2_ENGINE_GET_TYPE_F2C)(char *type, adios2_engine **engine, int *ierr)
{
    size_t size;
    *ierr = static_cast<int>(adios2_engine_get_type(type, &size, *engine));
}

#ifdef __cplusplus
}
#endif
