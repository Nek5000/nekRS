/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP3Base.h  base class for BP3Serializer and BP3Deserializer
 *
 *  Created on: Feb 2, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BP3_BP3BASE_H_
#define ADIOS2_TOOLKIT_FORMAT_BP_BP3_BP3BASE_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <bitset>
#include <map>
#include <memory> //std::shared_ptr
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/Engine.h"
#include "adios2/core/VariableBase.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/toolkit/aggregator/mpi/MPIChain.h"
#include "adios2/toolkit/format/buffer/heap/BufferSTL.h"
#include "adios2/toolkit/profiling/iochrono/IOChrono.h"

#include "adios2/toolkit/format/bp/BPBase.h"

namespace adios2
{
namespace format
{

/**
 * Base class for BP1Writer and BP1Reader format
 */
class BP3Base : virtual public BPBase
{

public:
    /**
     * Unique constructor
     * @param comm for m_BP1Aggregator
     */
    BP3Base(helper::Comm const &comm);

    virtual ~BP3Base() = default;

    /**
     * Vector version of BPBaseNames
     * @param names
     * @return vector of base (name.bp) names
     */
    std::vector<std::string> GetBPBaseNames(const std::vector<std::string> &names) const noexcept;

    /**
     * Get BP substream names from base names:
     * /path/name.bp.dir/name.bp.Index
     * where Index = Rank, and in aggregation = SubStreamID
     * @param baseNames inputs
     * @return vector of BP substream names for transports
     */
    std::vector<std::string>
    GetBPSubStreamNames(const std::vector<std::string> &baseNames) const noexcept;

    std::vector<std::string>
    GetBPMetadataFileNames(const std::vector<std::string> &names) const noexcept;

    std::string GetBPMetadataFileName(const std::string &name) const noexcept;

    std::string GetBPSubFileName(const std::string &name, const size_t subFileIndex,
                                 const bool hasSubFiles = true,
                                 const bool isReader = false) const noexcept;

    /**
     * Returns the estimated variable index size. Used by ResizeBuffer public
     * function
     * @param variableName input
     * @param count input variable local dimensions
     */
    size_t GetBPIndexSizeInData(const std::string &variableName, const Dims &count) const noexcept;

protected:
    ElementIndexHeader
    ReadElementIndexHeader(const std::vector<char> &buffer, size_t &position,
                           const bool isLittleEndian = true) const noexcept final;

private:
    std::string GetBPSubStreamName(const std::string &name, const size_t id,
                                   const bool hasSubFiles = true,
                                   const bool isReader = false) const noexcept;
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BP_BP3_BP3BASE_H_ */
