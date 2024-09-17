/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CampaignData.h
 * Campaign data from database
 *
 *  Created on: May 16, 2023
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#ifndef ADIOS2_ENGINE_CAMPAIGNDATA_H_
#define ADIOS2_ENGINE_CAMPAIGNDATA_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <sqlite3.h>

namespace adios2
{
namespace core
{
namespace engine
{

struct CampaignHost
{
    std::string hostname;
    std::string longhostname;
    std::vector<std::string> directory;
};

struct CampaignBPFile
{
    std::string name;
    size_t bpDatasetIdx; // index of parent CampaignBPDataset in the map
    bool compressed;
    size_t lengthOriginal;
    size_t lengthCompressed;
    int64_t ctime;
};

struct CampaignBPDataset
{
    std::string name;
    size_t hostIdx;
    size_t dirIdx;
    std::vector<CampaignBPFile> files;
};

struct CampaignData
{
    std::string version;
    std::vector<CampaignHost> hosts;
    std::map<size_t, CampaignBPDataset> bpdatasets;
};

void ReadCampaignData(sqlite3 *db, CampaignData &cd);

void SaveToFile(sqlite3 *db, const std::string &path, const CampaignBPFile &bpfile);

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_CAMPAIGDATA_H_ */
