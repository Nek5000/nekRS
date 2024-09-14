#ifndef ADIOS2_DERIVED_Function_TCC_
#define ADIOS2_DERIVED_Function_TCC_

#include "Function.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <numeric>

namespace adios2
{
namespace derived
{

template <class T>
T *ApplyOneToOne(std::vector<DerivedData> inputData, size_t dataSize,
                 std::function<T(T, T)> compFct)
{
    T *outValues = (T *)malloc(dataSize * sizeof(T));
    if (outValues == nullptr)
    {
        std::cout << "Allocation failed for the derived data" << std::endl;
        // TODO - throw an exception
    }
    memset((void *)outValues, 0, dataSize * sizeof(T));
    for (auto &variable : inputData)
    {
        for (size_t i = 0; i < dataSize; i++)
        {
            T data = *(reinterpret_cast<T *>(variable.Data) + i);
            outValues[i] = compFct(outValues[i], data);
        }
    }
    return outValues;
}

}
}
#endif
