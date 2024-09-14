/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_f2c_variable.cpp
 *
 *  Created on: Nov 12, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adios2_f2c_common.h"

#include "adios2/helper/adiosFunctions.h"

#ifdef __cplusplus
extern "C" {
#endif

void FC_GLOBAL(adios2_variable_name_f2c,
               ADIOS2_VARIABLE_NAME_F2C)(char *name, const adios2_variable **variable, int *ierr)
{
    size_t sizeC;
    *ierr = static_cast<int>(adios2_variable_name(name, &sizeC, *variable));
}

void FC_GLOBAL(adios2_variable_name_length_f2c,
               ADIOS2_VARIABLE_NAME_LENGTH_F2C)(int *size, const adios2_variable **variable,
                                                int *ierr)
{
    *size = -1;
    size_t sizeC;
    *ierr = static_cast<int>(adios2_variable_name(nullptr, &sizeC, *variable));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *size = static_cast<int>(sizeC);
    }
}

void FC_GLOBAL(adios2_variable_type_f2c,
               ADIOS2_VARIABLE_TYPE_F2C)(int *type, const adios2_variable **variable, int *ierr)
{
    *type = -1;
    adios2_type typeC;
    *ierr = static_cast<int>(adios2_variable_type(&typeC, *variable));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *type = static_cast<int>(typeC);
    }
}

void FC_GLOBAL(adios2_variable_ndims_f2c,
               ADIOS2_VARIABLE_NDIMS_F2C)(int *ndims, const adios2_variable **variable, int *ierr)
{
    *ndims = -1;
    size_t ndimsC;
    *ierr = static_cast<int>(adios2_variable_ndims(&ndimsC, *variable));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *ndims = static_cast<int>(ndimsC);
    }
}

void FC_GLOBAL(adios2_variable_shape_f2c,
               ADIOS2_VARIABLE_SHAPE_F2C)(int64_t *shape, const adios2_variable **variable,
                                          int *ierr)
{
    size_t ndims;
    *ierr = static_cast<int>(adios2_variable_ndims(&ndims, *variable));
    if (*ierr > 0)
    {
        return;
    }

    size_t shapeC[ndims];
    *ierr = static_cast<int>(adios2_variable_shape(shapeC, *variable));
    if (*ierr > 0)
    {
        return;
    }

    for (size_t d = 0; d < ndims; ++d)
    {
        shape[d] = static_cast<int64_t>(shapeC[d]);
    }
}

void FC_GLOBAL(adios2_variable_steps_f2c,
               adios2_variable_STEPS_F2C)(int64_t *steps, const adios2_variable **variable,
                                          int *ierr)
{
    *steps = -1;
    size_t stepsC;
    *ierr = static_cast<int>(adios2_variable_steps(&stepsC, *variable));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *steps = static_cast<int64_t>(stepsC);
    }
}

void FC_GLOBAL(adios2_set_memory_space_f2c, ADIOS2_SET_MEMORY_SPACE_F2C)(adios2_variable **variable,
                                                                         const int *mem, int *ierr)
{
    *ierr = static_cast<int>(
        adios2_set_memory_space(*variable, static_cast<adios2_memory_space>(*mem)));
}

void FC_GLOBAL(adios2_get_memory_space_f2c,
               ADIOS2_GET_MEMORY_SPACE_F2C)(int *mem, adios2_variable **variable, int *ierr)
{
    *mem = 0;
    adios2_memory_space Cmem;
    *ierr = static_cast<int>(adios2_get_memory_space(&Cmem, *variable));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *mem = static_cast<int>(Cmem);
    }
}

void FC_GLOBAL(adios2_set_shape_f2c, ADIOS2_SET_SHAPE_F2C)(adios2_variable **variable,
                                                           const int *ndims, const int64_t *shape,
                                                           int *ierr)
{
    auto lf_IntToSizeT = [](const int64_t *dimensions, const int size,
                            std::vector<std::size_t> &output) {
        output.resize(size);

        for (auto d = 0; d < size; ++d)
        {
            output[d] = dimensions[d];
        }
    };

    try
    {
        if (shape == nullptr || ndims == nullptr)
        {
            throw std::invalid_argument(
                "ERROR: either shape_dims, or ndims is a null pointer, in call "
                "to adios2_set_shape\n");
        }
        std::vector<std::size_t> shapeV;
        lf_IntToSizeT(shape, *ndims, shapeV);
        *ierr = static_cast<int>(adios2_set_shape(*variable, *ndims, shapeV.data()));
    }
    catch (...)
    {
        *ierr = static_cast<int>(adios2::helper::ExceptionToError("adios2_set_shape"));
    }
}

void FC_GLOBAL(adios2_set_block_selection_f2c,
               ADIOS2_SET_BLOCK_SELECTION_F2C)(adios2_variable **variable, const int64_t *block_id,
                                               int *ierr)
{
    *ierr = static_cast<int>(adios2_set_block_selection(*variable, *block_id));
}

void FC_GLOBAL(adios2_set_selection_f2c,
               ADIOS2_SET_SELECTION_F2C)(adios2_variable **variable, const int *ndims,
                                         const int64_t *start, const int64_t *count, int *ierr)
{
    auto lf_IntToSizeT = [](const int64_t *input, const int size) -> std::vector<std::size_t> {
        std::vector<std::size_t> output;
        output.reserve(size);

        for (auto d = 0; d < size; ++d)
        {
            output.push_back(static_cast<std::size_t>(input[d]));
        }
        return output;
    };

    try
    {
        // is this even possible coming from Fortran ?
        if (start == nullptr || count == nullptr || ndims == nullptr)
        {
            throw std::invalid_argument("ERROR: either start_dims, count_dims "
                                        "or ndims is a null pointer, in call "
                                        "to adios2_set_selection\n");
        }
        if (*ndims <= 0)
        {
            throw std::invalid_argument("ERROR: ndims must be larger than 0, "
                                        "in call to adios2_set_selection\n");
        }

        const std::vector<size_t> countV = lf_IntToSizeT(count, *ndims);
        // local array selection
        if (start[0] == -1)
        {
            *ierr =
                static_cast<int>(adios2_set_selection(*variable, *ndims, nullptr, countV.data()));
        }
        // global array selection
        else
        {
            const std::vector<std::size_t> startV = lf_IntToSizeT(start, *ndims);
            *ierr = static_cast<int>(
                adios2_set_selection(*variable, *ndims, startV.data(), countV.data()));
        }
    }
    catch (...)
    {
        *ierr = static_cast<int>(adios2::helper::ExceptionToError("adios2_set_selection"));
    }
}

void FC_GLOBAL(adios2_set_memory_selection_f2c,
               ADIOS2_SET_MEMORY_SELECTION_F2C)(adios2_variable **variable, const int *ndims,
                                                const int64_t *memory_start,
                                                const int64_t *memory_count, int *ierr)
{
    auto lf_IntToSizeT = [](const int64_t *dimensions, const int size,
                            std::vector<std::size_t> &output) {
        output.resize(size);

        for (auto d = 0; d < size; ++d)
        {
            output[d] = dimensions[d];
        }
    };

    try
    {
        if (memory_start == nullptr || memory_count == nullptr || ndims == nullptr)
        {
            throw std::invalid_argument("ERROR: either start_dims, count_dims "
                                        "or ndims is a null pointer, in call "
                                        "to adios2_set_memory_selection\n");
        }
        std::vector<std::size_t> memoryStartV, memoryCountV;
        lf_IntToSizeT(memory_start, *ndims, memoryStartV);
        lf_IntToSizeT(memory_count, *ndims, memoryCountV);
        *ierr = static_cast<int>(adios2_set_memory_selection(*variable, *ndims, memoryStartV.data(),
                                                             memoryCountV.data()));
    }
    catch (...)
    {
        *ierr = static_cast<int>(adios2::helper::ExceptionToError("adios2_set_memory_selection"));
    }
}

void FC_GLOBAL(adios2_set_step_selection_f2c,
               ADIOS2_SET_STEP_SELECTION_F2C)(adios2_variable **variable, const int64_t *step_start,
                                              const int64_t *step_count, int *ierr)
{
    try
    {
        if (step_start == nullptr || step_count == nullptr)
        {
            throw std::invalid_argument("ERROR: either step_start or step_count "
                                        "are null pointers, in call to "
                                        "adios2_set_step_selection\n");
        }

        if (step_start[0] < 0)
        {
            throw std::invalid_argument("ERROR: negative step_start in call to "
                                        "adios2_set_step_selection\n");
        }

        if (step_count[0] < 0)
        {
            throw std::invalid_argument("ERROR: negative step_count in call to "
                                        "adios2_set_step_selection\n");
        }

        const std::size_t stepStart = static_cast<std::size_t>(*step_start);
        const std::size_t stepCount = static_cast<std::size_t>(*step_count);
        *ierr = adios2_set_step_selection(*variable, stepStart, stepCount);
    }
    catch (...)
    {
        *ierr = static_cast<int>(adios2::helper::ExceptionToError("adios2_set_selection"));
    }
}

void FC_GLOBAL(adios2_add_operation_f2c,
               ADIOS2_ADD_OPERATION_F2C)(int *operation_id, adios2_variable **variable,
                                         adios2_operator **op, const char *key, const char *value,
                                         int *ierr)
{
    *operation_id = -1;
    size_t operation_idC;
    *ierr = static_cast<int>(adios2_add_operation(&operation_idC, *variable, *op, key, value));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *operation_id = static_cast<int>(operation_idC);
    }
}

void FC_GLOBAL(adios2_set_operation_parameter_f2c,
               ADIOS2_SET_OPERATION_PARAMETER_F2C)(adios2_variable **variable,
                                                   const int *operation_id, const char *key,
                                                   const char *value, int *ierr)
{
    try
    {
        if (*operation_id < 0)
        {
            throw std::invalid_argument("ERROR: operation_id can't be "
                                        "negative, in call to "
                                        "adios2_set_operation_paramter");
        }

        *ierr = static_cast<int>(adios2_set_operation_parameter(
            *variable, static_cast<std::size_t>(*operation_id), key, value));
    }
    catch (...)
    {
        *ierr =
            static_cast<int>(adios2::helper::ExceptionToError("adios2_set_operation_parameter"));
    }
}

void FC_GLOBAL(adios2_remove_operations_f2c,
               ADIOS2_REMOVE_OPERATIONS_F2C)(adios2_variable **variable, int *ierr)
{
    *ierr = static_cast<int>(adios2_remove_operations(*variable));
}

void FC_GLOBAL(adios2_variable_min_f2c,
               ADIOS2_VARIABLE_MIN_F2C)(void *min, const adios2_variable **variable, int *ierr)
{
    *ierr = static_cast<int>(adios2_variable_min(min, *variable));
}

void FC_GLOBAL(adios2_variable_max_f2c,
               ADIOS2_VARIABLE_MAX_F2C)(void *max, const adios2_variable **variable, int *ierr)
{
    *ierr = static_cast<int>(adios2_variable_max(max, *variable));
}

#ifdef __cplusplus
}
#endif
