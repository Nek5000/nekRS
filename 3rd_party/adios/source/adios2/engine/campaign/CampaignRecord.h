/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CampaignRecord.h
 * Campaign record struct
 *
 *  Created on: May 16, 2023
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#ifndef ADIOS2_ENGINE_CAMPAIGNRECORD_H_
#define ADIOS2_ENGINE_CAMPAIGNRECORD_H_

#include <string>
#include <unordered_map>
#include <vector>

namespace adios2
{
namespace core
{
namespace engine
{

/*
 * Data from processes to be recorded
 */

struct CampaignRecord
{
    std::vector<size_t> steps;
    size_t delta_step;
    std::vector<double> times;
    double delta_time;
    bool varying_deltas;

    CampaignRecord() : delta_step(0), delta_time(0.0), varying_deltas(false){};

    CampaignRecord(size_t step, double time) : delta_step(0), delta_time(0.0), varying_deltas(false)
    {
        steps.push_back(step);
        times.push_back(time);
    };
};

using CampaignRecordMap = std::unordered_map<std::string, CampaignRecord>;

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_CAMPAIGRECORD_H_ */
