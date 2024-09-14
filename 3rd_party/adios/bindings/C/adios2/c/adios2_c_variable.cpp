/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_variable.cpp
 *
 *  Created on: Nov 10, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adios2_c_variable.h"

#include "adios2/core/Variable.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2_c_internal.h"

namespace
{
adios2_shapeid adios2_ToShapeID(const adios2::ShapeID shapeIDCpp, const std::string &hint)
{
    adios2_shapeid shapeID = adios2_shapeid_unknown;
    switch (shapeIDCpp)
    {

    case (adios2::ShapeID::GlobalValue):
        shapeID = adios2_shapeid_global_value;
        break;

    case (adios2::ShapeID::GlobalArray):
        shapeID = adios2_shapeid_global_array;
        break;

    case (adios2::ShapeID::JoinedArray):
        shapeID = adios2_shapeid_joined_array;
        break;

    case (adios2::ShapeID::LocalValue):
        shapeID = adios2_shapeid_local_value;
        break;

    case (adios2::ShapeID::LocalArray):
        shapeID = adios2_shapeid_local_array;
        break;

    default:
        throw std::invalid_argument("ERROR: invalid adios2_shapeid, " + hint + "\n");
    }

    return shapeID;
}

} // end anonymous namespace

#ifdef __cplusplus
extern "C" {
#endif

adios2_error adios2_set_shape(adios2_variable *variable, const size_t ndims, const size_t *shape)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call to "
                                                  "adios2_set_shape");
        adios2::helper::CheckForNullptr(shape, "for start, in call to "
                                               "adios2_set_shape");

        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);

        const adios2::Dims shapeV(shape, shape + ndims);
        variableBase->SetShape(shapeV);
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_set_shape"));
    }
}

adios2::MemorySpace adios2_ToMemorySpace(const adios2_memory_space Cmem)
{
    adios2::MemorySpace mem = adios2::MemorySpace::Detect;
    switch (Cmem)
    {

    case adios2_memory_space_host:
        mem = adios2::MemorySpace::Host;
        break;

#ifdef ADIOS2_HAVE_GPU_SUPPORT
    case adios2_memory_space_gpu:
        mem = adios2::MemorySpace::GPU;
        break;
#endif
    default:
        break;
    }
    return mem;
}

adios2_memory_space adios2_FromMemorySpace(const adios2::MemorySpace mem)
{
    adios2_memory_space Cmem = adios2_memory_space_detect;
    switch (mem)
    {

    case adios2::MemorySpace::Host:
        Cmem = adios2_memory_space_host;
        break;

#ifdef ADIOS2_HAVE_GPU_SUPPORT
    case adios2::MemorySpace::GPU:
        Cmem = adios2_memory_space_gpu;
        break;
#endif
    default:
        break;
    }
    return Cmem;
}

adios2_error adios2_set_memory_space(adios2_variable *variable, const adios2_memory_space mem)
{
    try
    {
        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);
        variableBase->SetMemorySpace(adios2_ToMemorySpace(mem));

        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(
            adios2::helper::ExceptionToError("adios2_set_memory_space"));
    }
}

adios2_error adios2_get_memory_space(adios2_memory_space *mem, adios2_variable *variable)
{
    try
    {
        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);
        *mem = adios2_FromMemorySpace(variableBase->m_MemSpace);

        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(
            adios2::helper::ExceptionToError("adios2_set_memory_space"));
    }
}

adios2_error adios2_set_block_selection(adios2_variable *variable, const size_t block_id)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call to "
                                                  "adios2_set_block_selection");

        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);
        variableBase->SetBlockSelection(block_id);
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(
            adios2::helper::ExceptionToError("adios2_set_block_selection"));
    }
}

adios2_error adios2_set_selection(adios2_variable *variable, const size_t ndims,
                                  const size_t *start, const size_t *count)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call to "
                                                  "adios2_set_selection");
        adios2::helper::CheckForNullptr(count, "for count, in call to "
                                               "adios2_set_selection");
        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);

        const adios2::Dims startV =
            (start == nullptr) ? adios2::Dims() : adios2::Dims(start, start + ndims);
        const adios2::Dims countV = adios2::Dims(count, count + ndims);

        variableBase->SetSelection({startV, countV});
        variableBase->CheckDimensions("in call to adios2_set_selection");
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_set_selection"));
    }
}

adios2_error adios2_set_memory_selection(adios2_variable *variable, const size_t ndims,
                                         const size_t *memory_start, const size_t *memory_count)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call to "
                                                  "adios2_set_memory_selection");
        adios2::helper::CheckForNullptr(memory_start, "for start, in call to "
                                                      "adios2_set_memory_selection");
        adios2::helper::CheckForNullptr(memory_count, "for count, in call to "
                                                      "adios2_set_memory_selection");
        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);

        const adios2::Dims memoryStartV(memory_start, memory_start + ndims);
        const adios2::Dims memoryCountV(memory_count, memory_count + ndims);
        variableBase->SetMemorySelection({memoryStartV, memoryCountV});
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(
            adios2::helper::ExceptionToError("adios2_set_memory_selection"));
    }
}

adios2_error adios2_set_step_selection(adios2_variable *variable, const size_t step_start,
                                       const size_t step_count)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call to "
                                                  "adios2_set_step_selection");
        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);
        variableBase->SetStepSelection(adios2::Box<size_t>{step_start, step_count});
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_set_selection"));
    }
}

adios2_error adios2_variable_name(char *name, size_t *size, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(
            variable, "for const adios2_variable, in call to adios2_variable_name");
        adios2::helper::CheckForNullptr(size, "for size_t* size, in call to adios2_variable_name");

        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);
        return String2CAPI(variableBase->m_Name, name, size);
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_variable_name"));
    }
}

adios2_error adios2_variable_type(adios2_type *type, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(
            variable, "for const adios2_variable, in call to adios2_variable_type");

        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);

        const adios2::DataType typeCpp = variableBase->m_Type;
        if (typeCpp == adios2::helper::GetDataType<std::string>())
        {
            *type = adios2_type_string;
        }
#define make_case(T)                                                                               \
    else if (typeCpp == adios2::helper::GetDataType<MapAdios2Type<T>::Type>())                     \
    {                                                                                              \
        *type = T;                                                                                 \
    }
        ADIOS2_FOREACH_C_TYPE_1ARG(make_case)
#undef make_case
        else
        {
            *type = adios2_type_unknown;
        }
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_variable_type"));
    }
}

adios2_error adios2_variable_type_string(char *type, size_t *size, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for const adios2_variable, in call to "
                                                  "adios2_variable_type_string");
        adios2::helper::CheckForNullptr(
            size, "for size_t* length, in call to adios2_variable_type_string");

        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);
        return String2CAPI(ToString(variableBase->m_Type), type, size);
    }
    catch (...)
    {
        return static_cast<adios2_error>(
            adios2::helper::ExceptionToError("adios2_variable_type_string"));
    }
}

adios2_error adios2_variable_shapeid(adios2_shapeid *shapeid, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for const adios2_variable, in call to "
                                                  "adios2_variable_shapeid");
        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);

        *shapeid = adios2_ToShapeID(variableBase->m_ShapeID, "in call to adios2_variable_shapeid");
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(
            adios2::helper::ExceptionToError("adios2_variable_shapeid"));
    }
}

adios2_error adios2_variable_ndims(size_t *ndims, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for const adios2_variable, in call to "
                                                  "adios2_variable_ndims");
        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);
        *ndims = variableBase->m_Count.size();
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_variable_ndims"));
    }
}

adios2_error adios2_variable_shape_with_memory_space(size_t *shape, const adios2_variable *variable,
                                                     adios2_memory_space mem)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for const adios2_variable, in call to "
                                                  "adios2_variable_shape");
        adios2::helper::CheckForNullptr(shape, "for size_t* shape, in call to "
                                               "adios2_variable_shape");

        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);

        const adios2::DataType typeCpp = variableBase->m_Type;
        if (typeCpp == adios2::DataType::Struct)
        {
            // not supported
        }
#define declare_template_instantiation(T)                                                          \
    else if (typeCpp == adios2::helper::GetDataType<T>())                                          \
    {                                                                                              \
        const adios2::core::Variable<T> *variable =                                                \
            dynamic_cast<const adios2::core::Variable<T> *>(variableBase);                         \
        adios2::Dims shapeCpp;                                                                     \
        if (mem == adios2_memory_space_host)                                                       \
            shapeCpp = variable->Shape(adios2::EngineCurrentStep);                                 \
        else                                                                                       \
            shapeCpp = variable->Shape(adios2::EngineCurrentStep, adios2_ToMemorySpace(mem));      \
        std::copy(shapeCpp.begin(), shapeCpp.end(), shape);                                        \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_variable_shape"));
    }
}

adios2_error adios2_variable_shape(size_t *shape, const adios2_variable *variable)
{
    return adios2_variable_shape_with_memory_space(shape, variable, adios2_memory_space_host);
}

adios2_error adios2_variable_start(size_t *start, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for const adios2_variable, in call to "
                                                  "adios2_variable_start");
        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);

        std::copy(variableBase->m_Start.begin(), variableBase->m_Start.end(), start);
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_variable_shape"));
    }
}

adios2_error adios2_variable_count(size_t *count, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for const adios2_variable, in call to "
                                                  "adios2_variable_count");
        adios2::helper::CheckForNullptr(count, "for const adios2_count, in call to "
                                               "adios2_variable_count");

        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);

        const adios2::DataType typeCpp = variableBase->m_Type;
        if (typeCpp == adios2::DataType::Struct)
        {
            // not supported
        }
#define declare_template_instantiation(T)                                                          \
    else if (typeCpp == adios2::helper::GetDataType<T>())                                          \
    {                                                                                              \
        const adios2::core::Variable<T> *variable =                                                \
            dynamic_cast<const adios2::core::Variable<T> *>(variableBase);                         \
        const adios2::Dims countCpp = variable->Count();                                           \
        std::copy(countCpp.begin(), countCpp.end(), count);                                        \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_variable_count"));
    }
}

adios2_error adios2_variable_steps_start(size_t *steps_start, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for const adios2_variable, in call to "
                                                  "adios2_variable_steps_start");
        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);
        *steps_start = variableBase->m_AvailableStepsStart;
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(
            adios2::helper::ExceptionToError("adios2_variable_steps_start"));
    }
}

adios2_error adios2_variable_steps(size_t *steps, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for const adios2_variable, in call to "
                                                  "adios2_variable_steps");
        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);
        *steps = variableBase->m_AvailableStepsCount;
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_variable_steps"));
    }
}

adios2_error adios2_selection_size(size_t *size, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call to "
                                                  "adios2_selection_size");
        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);

        const adios2::DataType typeCpp = variableBase->m_Type;

        if (typeCpp == adios2::DataType::Struct)
        {
            // not supported
        }
#define declare_template_instantiation(T)                                                          \
    else if (typeCpp == adios2::helper::GetDataType<T>())                                          \
    {                                                                                              \
        const adios2::core::Variable<T> *variable =                                                \
            dynamic_cast<const adios2::core::Variable<T> *>(variableBase);                         \
        *size = variable->SelectionSize();                                                         \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_selection_size"));
    }
}

adios2_error adios2_add_operation(size_t *operation_index, adios2_variable *variable,
                                  adios2_operator *op, const char *key, const char *value)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call to "
                                                  "adios2_add_operation");
        adios2::helper::CheckForNullptr(op, "for adios2_operator, in call to "
                                            "adios2_add_operation");

        adios2::helper::CheckForNullptr(key, "for char* key, in call to "
                                             "adios2_add_operation");

        adios2::helper::CheckForNullptr(value, "for char* value, in call to "
                                               "adios2_add_operation");

        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);

        auto *opCpp = reinterpret_cast<std::pair<std::string, adios2::Params> *>(op);

        auto params = adios2::Params{{key, value}};

        for (const auto &p : opCpp->second)
        {
            params[p.first] = p.second;
        }

        *operation_index = variableBase->AddOperation(opCpp->first, params);

        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_add_operation"));
    }
}

adios2_error adios2_set_operation_parameter(adios2_variable *variable, const size_t operation_id,
                                            const char *key, const char *value)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call to "
                                                  "adios2_set_operation_parameter");

        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);
        variableBase->SetOperationParameter(operation_id, key, value);
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(
            adios2::helper::ExceptionToError("adios2_set_operation_parameter"));
    }
}

adios2_error adios2_remove_operations(adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call to "
                                                  "adios2_remove_operations");
        adios2::core::VariableBase *variableBase =
            reinterpret_cast<adios2::core::VariableBase *>(variable);
        variableBase->RemoveOperations();
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(
            adios2::helper::ExceptionToError("adios2_remove_operations"));
    }
}

adios2_error adios2_variable_min(void *min, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call "
                                                  "to adios2_variable_min");
        adios2::helper::CheckForNullptr(min, "for void* min, in call to adios2_variable_min");

        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);
        const adios2::DataType type(variableBase->m_Type);

        if (type == adios2::DataType::Struct)
        {
            // not supported
        }
#define declare_template_instantiation(T)                                                          \
    else if (type == adios2::helper::GetDataType<T>())                                             \
    {                                                                                              \
        T *minT = reinterpret_cast<T *>(min);                                                      \
        const adios2::core::Variable<T> *variableT =                                               \
            dynamic_cast<const adios2::core::Variable<T> *>(variableBase);                         \
        *minT = variableT->Min(adios2::EngineCurrentStep);                                         \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_variable_min"));
    }
}

adios2_error adios2_variable_max(void *max, const adios2_variable *variable)
{
    try
    {
        adios2::helper::CheckForNullptr(variable, "for adios2_variable, in call "
                                                  "to adios2_variable_max");
        adios2::helper::CheckForNullptr(max, "for void* max, in call to adios2_variable_max");

        const adios2::core::VariableBase *variableBase =
            reinterpret_cast<const adios2::core::VariableBase *>(variable);
        const adios2::DataType type(variableBase->m_Type);

        if (type == adios2::DataType::Struct)
        {
            // not supported
        }
#define declare_template_instantiation(T)                                                          \
    else if (type == adios2::helper::GetDataType<T>())                                             \
    {                                                                                              \
        T *maxT = reinterpret_cast<T *>(max);                                                      \
        const adios2::core::Variable<T> *variableT =                                               \
            dynamic_cast<const adios2::core::Variable<T> *>(variableBase);                         \
        *maxT = variableT->Max(adios2::EngineCurrentStep);                                         \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
        return adios2_error_none;
    }
    catch (...)
    {
        return static_cast<adios2_error>(adios2::helper::ExceptionToError("adios2_variable_max"));
    }
}

#ifdef __cplusplus
} // end extern C
#endif
