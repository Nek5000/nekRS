/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPOperation.h :
 *
 *  Created on: Jul 12, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPOPERATION_BPOPERATION_H_
#define ADIOS2_TOOLKIT_FORMAT_BP_BPOPERATION_BPOPERATION_H_

#include <string>
#include <vector>

#include "adios2/common/ADIOSMacros.h"
#include "adios2/core/Variable.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/toolkit/format/buffer/heap/BufferSTL.h"

namespace adios2
{
namespace format
{

class BPBackCompatOperation
{
public:
    BPBackCompatOperation() = default;
    virtual ~BPBackCompatOperation() = default;

    /**
     * Deserializes metadata in the form of parameters
     * @param buffer contains serialized metadata buffer
     * @param info parameters info from metadata buffer
     */
    virtual void GetMetadata(const std::vector<char> &buffer, Params &info) const noexcept = 0;

    virtual void GetData(const char *input, const helper::BlockOperationInfo &blockOperationInfo,
                         char *dataOutput) const = 0;
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BP_OPERATION_BPOPERATION_H_ */
