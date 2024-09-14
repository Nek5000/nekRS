/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * OperatorFactory.h :
 *
 *  Created on: Sep 29, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/Operator.h"
#include <memory>

namespace adios2
{
namespace core
{

std::string OperatorTypeToString(const Operator::OperatorType type);

std::shared_ptr<Operator> MakeOperator(const std::string &type, const Params &parameters);

size_t Decompress(const char *bufferIn, const size_t sizeIn, char *dataOut, MemorySpace memSpace,
                  std::shared_ptr<Operator> op = nullptr);

} // end namespace core
} // end namespace adios2
