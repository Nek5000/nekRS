/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_f2c_io.cpp
 *
 *  Created on: Nov 8, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adios2_f2c_common.h"

#include <cstdint>
#include <stdexcept> //std::invalid_argument

#include <string.h> //strcpy

#include "adios2/common/ADIOSTypes.h"
#include "adios2/helper/adiosFunctions.h"

#ifdef __cplusplus
extern "C" {
#endif

void FC_GLOBAL(adios2_in_config_file_f2c,
               ADIOS2_IN_CONFIG_FILE_F2C)(int *result, const adios2_io **io, int *ierr)
{
    adios2_bool resultC;
    *ierr = static_cast<int>(adios2_in_config_file(&resultC, *io));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *result = (resultC == adios2_true) ? 1 : 0;
    }
}

void FC_GLOBAL(adios2_set_engine_f2c, ADIOS2_SET_ENGINE_F2C)(adios2_io **io,
                                                             const char *engine_type, int *ierr)
{
    *ierr = static_cast<int>(adios2_set_engine(*io, engine_type));
}

void FC_GLOBAL(adios2_set_parameters_f2c,
               ADIOS2_SET_PARAMETERS_F2C)(adios2_io **io, const char *parameters, int *ierr)
{
    *ierr = static_cast<int>(adios2_set_parameters(*io, parameters));
}

void FC_GLOBAL(adios2_set_parameter_f2c, ADIOS2_SET_PARAMETER_F2C)(adios2_io **io, const char *key,
                                                                   const char *value, int *ierr)
{
    *ierr = static_cast<int>(adios2_set_parameter(*io, key, value));
}

void FC_GLOBAL(adios2_get_parameter_f2c, ADIOS2_GET_PARAMETER_F2C)(char *value,
                                                                   const adios2_io **io,
                                                                   const char *key, int *ierr)
{
    size_t sizeC;
    *ierr = static_cast<int>(adios2_get_parameter(value, &sizeC, *io, key));
}

void FC_GLOBAL(adios2_get_parameter_length_f2c,
               ADIOS2_GET_PARAMETER_LENGTH_F2C)(int *size, const adios2_io **io, const char *key,
                                                int *ierr)
{
    *size = -1;
    size_t sizeC;
    *ierr = static_cast<int>(adios2_get_parameter(nullptr, &sizeC, *io, key));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *size = static_cast<int>(sizeC);
    }
}

void FC_GLOBAL(adios2_clear_parameters_f2c, ADIOS2_CLEAR_PARAMETERS_F2C)(adios2_io **io, int *ierr)
{
    *ierr = static_cast<int>(adios2_clear_parameters(*io));
}

void FC_GLOBAL(adios2_add_transport_f2c,
               ADIOS2_ADD_TRANSPORT_F2C)(int *transport_index, adios2_io **io,
                                         const char *transport_type, int *ierr)
{

    *transport_index = -1;
    std::size_t transport_indexC;
    *ierr = static_cast<int>(adios2_add_transport(&transport_indexC, *io, transport_type));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *transport_index = static_cast<int>(transport_indexC);
    }
}

void FC_GLOBAL(adios2_set_transport_parameter_f2c,
               ADIOS2_SET_TRANSPORT_PARAMETER_F2C)(adios2_io **io, const int *transport_index,
                                                   const char *key, const char *value, int *ierr)
{
    if (*transport_index < 0)
    {
        throw std::invalid_argument("ERROR: transport_index can't be negative, "
                                    "in call to "
                                    "adios2_set_transport_parameter");
    }
    const std::size_t transport_indexC = static_cast<std::size_t>(*transport_index);
    *ierr = static_cast<int>(adios2_set_transport_parameter(*io, transport_indexC, key, value));
}

// needed due to lack of nullptr in Fortran
void FC_GLOBAL(adios2_define_global_variable_f2c,
               ADIOS2_DEFINE_GLOBAL_VARIABLE_F2C)(adios2_variable **variable, adios2_io **io,
                                                  const char *name, const int *type, int *ierr)
{
    *variable = adios2_define_variable(*io, name, static_cast<adios2_type>(*type), 0, nullptr,
                                       nullptr, nullptr, adios2_constant_dims_true);
    *ierr = (*variable == NULL) ? static_cast<int>(adios2_error_exception)
                                : static_cast<int>(adios2_error_none);
}

void FC_GLOBAL(adios2_define_variable_f2c,
               ADIOS2_DEFINE_VARIABLE_F2C)(adios2_variable **variable, adios2_io **io,
                                           const char *name, const int *type, const int *ndims,
                                           const int64_t *shape, const int64_t *start,
                                           const int64_t *count, const int *constant_dims,
                                           int *ierr)
{
    auto lf_IntToSizeT = [](const int64_t *dimensions, const int size,
                            const std::string hint) -> adios2::Dims {
        adios2::Dims output(static_cast<size_t>(size));

        for (auto d = 0; d < size; ++d)
        {
            if (dimensions[d] < 0)
            {
                throw std::invalid_argument("ERROR: negative dimension " + hint +
                                            ", in call to adios2_define_variable");
            }
            output[d] = dimensions[d];
        }
        return output;
    };

    try
    {
        if (*ndims <= 0)
        {
            throw std::invalid_argument("ERROR: negative ndims in Fortran, in call to "
                                        "adios2_define_variable");
        }

        // Check for local value
        if (*ndims == 1)
        {
            if (shape[0] == -2)
            {
                size_t shapeC[1];
                shapeC[0] = adios2_local_value_dim;

                *variable = adios2_define_variable(
                    *io, name, static_cast<adios2_type>(*type), *ndims, shapeC, nullptr, nullptr,
                    static_cast<adios2_constant_dims>(*constant_dims));
                *ierr = (*variable == NULL) ? static_cast<int>(adios2_error_exception)
                                            : static_cast<int>(adios2_error_none);
                return;
            }
        }
        // Check for local variables
        if (shape[0] == -1)
        {
            if (start[0] != -1)
            {
                throw std::invalid_argument("ERROR: both shape and start must be "
                                            "adios2_null_dims when declaring "
                                            "local variables in Fortran, in call to "
                                            "adios2_define_variable");
            }

            const adios2::Dims countV = lf_IntToSizeT(count, *ndims, "count");

            *variable = adios2_define_variable(*io, name, static_cast<adios2_type>(*type), *ndims,
                                               nullptr, nullptr, countV.data(),
                                               static_cast<adios2_constant_dims>(*constant_dims));
            *ierr = (*variable == NULL) ? static_cast<int>(adios2_error_exception)
                                        : static_cast<int>(adios2_error_none);
            return;
        }

        const adios2::Dims shapeV = lf_IntToSizeT(shape, *ndims, "shape");
        const adios2::Dims startV = lf_IntToSizeT(start, *ndims, "start");
        const adios2::Dims countV = lf_IntToSizeT(count, *ndims, "count");

        *variable = adios2_define_variable(*io, name, static_cast<adios2_type>(*type),
                                           static_cast<std::size_t>(*ndims), shapeV.data(),
                                           startV.data(), countV.data(),
                                           static_cast<adios2_constant_dims>(*constant_dims));
        *ierr = (*variable == NULL) ? static_cast<int>(adios2_error_exception)
                                    : static_cast<int>(adios2_error_none);
    }
    catch (std::exception &e)
    {
        *ierr = static_cast<int>(adios2::helper::ExceptionToError("adios2_define_variable"));
    }
}

void FC_GLOBAL(adios2_define_derived_variable_f2c,
               ADIOS2_DEFINE_DERIVED_VARIABLE_F2C)(adios2_derived_variable **variable,
                                                   adios2_io **io, const char *name,
                                                   const char *expression, const int *type,
                                                   int *ierr)
{
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
    *variable = adios2_define_derived_variable(*io, name, expression,
                                               static_cast<adios2_derived_var_type>(*type));
    *ierr = (*variable == NULL) ? static_cast<int>(adios2_error_exception)
                                : static_cast<int>(adios2_error_none);
#else
    std::cout << "ADIOS2 Warning:  adios2_define_derived_variable() is not supported in the "
                 "current ADIOS2 build. The expression "
              << expression << " will be ignored and the variable " << name
              << " will not be produced." << std::endl;
    *variable = nullptr;
    *ierr = static_cast<int>(adios2_error_exception);
#endif
}

struct cnamelist
{
    char **names;
    size_t count;
};

void FC_GLOBAL(adios2_available_variables_f2c,
               ADIOS2_AVAILABLE_VARIABLES_F2C)(adios2_io **io, int64_t *namestruct, int *vars_count,
                                               int *max_var_name_len, int *ierr)
{
    cnamelist *info = new (cnamelist);
    info->names = adios2_available_variables(*io, &info->count);
    *vars_count = static_cast<int>(info->count);
    size_t maxlen = 0;
    for (size_t i = 0; i < info->count; ++i)
    {
        auto l = strlen(info->names[i]);
        if (l > maxlen)
        {
            maxlen = l;
        }
    }
    *max_var_name_len = static_cast<int>(maxlen);
    *namestruct = static_cast<int64_t>(reinterpret_cast<std::uintptr_t>(info));
    *ierr = 0;
}

void FC_GLOBAL(adios2_available_attributes_f2c,
               ADIOS2_AVAILABLE_ATTRIBUTES_F2C)(adios2_io **io, int64_t *namestruct,
                                                int *attrs_count, int *max_attr_name_len, int *ierr)
{
    cnamelist *info = new (cnamelist);
    info->names = adios2_available_attributes(*io, &info->count);
    *attrs_count = static_cast<int>(info->count);

    size_t maxlen = 0;
    for (size_t i = 0; i < info->count; ++i)
    {
        auto l = strlen(info->names[i]);
        if (l > maxlen)
        {
            maxlen = l;
        }
    }
    *max_attr_name_len = static_cast<int>(maxlen);

    *namestruct = static_cast<int64_t>(reinterpret_cast<std::uintptr_t>(info));
    *ierr = 0;
}

void FC_GLOBAL(adios2_retrieve_namelist_f2c,
               ADIOS2_RETRIEVE_NAMELIST_F2C)(int64_t *namestruct, void *namelist, int *ierr,
                                             int namelist_len)
{
    cnamelist *info = reinterpret_cast<cnamelist *>(*namestruct);
    if (info != NULL)
    {
        for (size_t i = 0; i < info->count; i++)
        {
            char *fs = (char *)namelist + i * namelist_len;
            size_t len = strlen(info->names[i]);
            if (len > static_cast<size_t>(namelist_len))
            {
                len = static_cast<size_t>(namelist_len);
            }
            // copy C string without '\0'
            memcpy(fs, info->names[i], len);
            // pad with spaces
            memset(fs + len, ' ', namelist_len - len);
        }
        *ierr = 0;
        delete (info);
        *namestruct = 0;
    }
    else
    {
        *ierr = 1;
    }
}

void FC_GLOBAL(adios2_inquire_variable_f2c,
               ADIOS2_INQUIRE_VARIABLE_F2C)(adios2_variable **variable, adios2_io **io,
                                            const char *name, int *ierr)
{

    *variable = adios2_inquire_variable(*io, name);
    *ierr = (*variable == NULL) ? static_cast<int>(adios2_error_exception)
                                : static_cast<int>(adios2_error_none);
}

void FC_GLOBAL(adios2_remove_variable_f2c, ADIOS2_REMOVE_VARIABLE_F2C)(int *result, adios2_io **io,
                                                                       const char *name, int *ierr)
{
    adios2_bool resultC;
    *ierr = static_cast<int>(adios2_remove_variable(&resultC, *io, name));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *result = (resultC == adios2_true) ? 1 : 0;
    }
}

void FC_GLOBAL(adios2_remove_all_variables_f2c, ADIOS2_REMOVE_ALL_VARIABLES_F2C)(adios2_io **io,
                                                                                 int *ierr)
{
    *ierr = static_cast<int>(adios2_remove_all_variables(*io));
}

void FC_GLOBAL(adios2_define_vattr_f2c,
               ADIOS2_DEFINE_VATTR_F2C)(adios2_attribute **attribute, adios2_io **io,
                                        const char *name, const int *type, const void *value,
                                        const char *variable_name, const char *separator, int *ierr)
{

    *attribute = adios2_define_variable_attribute(*io, name, static_cast<adios2_type>(*type), value,
                                                  variable_name, separator);
    *ierr = (*attribute == NULL) ? static_cast<int>(adios2_error_exception)
                                 : static_cast<int>(adios2_error_none);
}

void FC_GLOBAL(adios2_define_vattr_array_f2c,
               ADIOS2_DEFINE_VATTR_ARRAY_F2C)(adios2_attribute **attribute, adios2_io **io,
                                              const char *name, const int *type, const void *data,
                                              const int *size, const char *variable_name,
                                              const char *separator, int *ierr)
{
    try
    {
        if (*size < 1)
        {
            throw std::invalid_argument("ERROR: attribute " + std::string(name) +
                                        " array size must be larger or equal to 1, in call to "
                                        "adios2_define_attribute\n");
        }

        if (*type == adios2_type_string)
        {
            char **char2D = new char *[*size];

            // need to covert to row-major char** style, from Fortran large
            // sequential array of character*(len=4096)
            for (auto i = 0; i < *size; ++i)
            {
                char2D[i] = new char[adios2_string_array_element_max_size];

                const char *fstringMemory = reinterpret_cast<const char *>(data);

                strcpy(char2D[i], &fstringMemory[i * adios2_string_array_element_max_size]);
            }

            *attribute = adios2_define_variable_attribute_array(
                *io, name, static_cast<adios2_type>(*type), static_cast<const void *>(char2D),
                static_cast<std::size_t>(*size), variable_name, separator);

            delete[] char2D;
        }
        else
        {
            *attribute = adios2_define_variable_attribute_array(
                *io, name, static_cast<adios2_type>(*type), data, static_cast<std::size_t>(*size),
                variable_name, separator);
        }
        *ierr = (*attribute == NULL) ? static_cast<int>(adios2_error_exception)
                                     : static_cast<int>(adios2_error_none);
    }
    catch (...)
    {
        *ierr = adios2_error_exception;
    }
}

void FC_GLOBAL(adios2_inquire_attribute_f2c,
               ADIOS2_INQUIRE_ATTRIBUTE_F2C)(adios2_attribute **attribute, adios2_io **io,
                                             const char *attribute_name, int *ierr)
{
    *attribute = adios2_inquire_attribute(*io, attribute_name);
    *ierr = (*attribute == NULL) ? static_cast<int>(adios2_error_exception)
                                 : static_cast<int>(adios2_error_none);
}

void FC_GLOBAL(adios2_inquire_variable_attribute_f2c,
               ADIOS2_INQUIRE_VARIABLE_ATTRIBUTE_F2C)(adios2_attribute **attribute, adios2_io **io,
                                                      const char *attribute_name,
                                                      const char *variable_name,
                                                      const char *separator, int *ierr)
{
    *attribute = adios2_inquire_variable_attribute(*io, attribute_name, variable_name, separator);
    *ierr = (*attribute == NULL) ? static_cast<int>(adios2_error_exception)
                                 : static_cast<int>(adios2_error_none);
}

void FC_GLOBAL(adios2_remove_attribute_f2c,
               ADIOS2_REMOVE_ATTRIBUTE_F2C)(int *result, adios2_io **io, const char *name,
                                            int *ierr)
{
    adios2_bool resultC;
    *ierr = static_cast<int>(adios2_remove_attribute(&resultC, *io, name));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *result = (resultC == adios2_true) ? 1 : 0;
    }
}

void FC_GLOBAL(adios2_remove_all_attributes_f2c, ADIOS2_REMOVE_ALL_ATTRIBUTES_F2C)(adios2_io **io,
                                                                                   int *ierr)
{
    *ierr = static_cast<int>(adios2_remove_all_attributes(*io));
}

void FC_GLOBAL(adios2_open_f2c, ADIOS2_OPEN_F2C)(adios2_engine **engine, adios2_io **io,
                                                 const char *name, const int *open_mode, int *ierr)
{
    *engine = adios2_open(*io, name, static_cast<adios2_mode>(*open_mode));
    *ierr = (*engine == NULL) ? static_cast<int>(adios2_error_exception)
                              : static_cast<int>(adios2_error_none);
}

void FC_GLOBAL(adios2_flush_all_engines_f2c, ADIOS2_FLUSH_ALL_ENGINES_F2C)(adios2_io **io,
                                                                           int *ierr)
{
    *ierr = static_cast<int>(adios2_flush_all_engines(*io));
}

void FC_GLOBAL(adios2_io_engine_type_f2c,
               ADIOS2_IO_ENGINE_TYPE_F2C)(char *type, const adios2_io **io, int *ierr)
{
    size_t sizeC;
    *ierr = static_cast<int>(adios2_engine_type(type, &sizeC, *io));
}

void FC_GLOBAL(adios2_io_engine_type_length_f2c,
               ADIOS2_io_ENGINE_TYPE_LENGTH_F2C)(int *size, const adios2_io **io, int *ierr)
{
    *size = -1;
    size_t sizeC;
    *ierr = static_cast<int>(adios2_engine_type(nullptr, &sizeC, *io));
    if (*ierr == static_cast<int>(adios2_error_none))
    {
        *size = static_cast<int>(sizeC);
    }
}

#ifdef __cplusplus
}
#endif
