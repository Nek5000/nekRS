#ifndef ADIOS2_DERIVED_Function_H_
#define ADIOS2_DERIVED_Function_H_

#include "ExprHelper.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/helper/adiosLog.h"
#include <functional>

namespace adios2
{
namespace derived
{

struct DerivedData
{
    void *Data;
    Dims Start;
    Dims Count;
};

struct OperatorFunctions
{
    std::function<DerivedData(std::vector<DerivedData>, DataType)> ComputeFct;
    std::function<Dims(std::vector<Dims>)> DimsFct;
};

DerivedData AddFunc(std::vector<DerivedData> input, DataType type);
DerivedData MagnitudeFunc(std::vector<DerivedData> input, DataType type);
DerivedData Curl3DFunc(std::vector<DerivedData> input, DataType type);

template <class T>
T linear_interp(T *data, size_t index, size_t count, size_t stride = 1);

Dims SameDimsFunc(std::vector<Dims> input);
Dims CurlDimsFunc(std::vector<Dims> input);

const std::map<adios2::detail::ExpressionOperator, OperatorFunctions> OpFunctions = {
    {adios2::detail::ExpressionOperator::OP_ADD, {AddFunc, SameDimsFunc}},
    {adios2::detail::ExpressionOperator::OP_CURL, {Curl3DFunc, CurlDimsFunc}},
    {adios2::detail::ExpressionOperator::OP_MAGN, {MagnitudeFunc, SameDimsFunc}}};

template <class T>
T *ApplyOneToOne(std::vector<DerivedData> inputData, size_t dataSize,
                 std::function<T(T, T)> compFct);

}
}
#endif
