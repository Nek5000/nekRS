/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CampaignManager.cpp
 *
 * This is NOT a writer Engine but the CampaignReader is a reader Engine.
 *
 *  Created on: May 15, 2023
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#include "CampaignManager.h"

#include "adios2/helper/adiosFunctions.h"

#include <iostream>

#include <nlohmann_json.hpp>
#include <sqlite3.h>

namespace adios2
{
namespace core
{
namespace engine
{

int CMapToSqlite(const CampaignRecordMap &cmap, const int rank, std::string name)
{
    sqlite3 *db;
    int rc;
    char *zErrMsg = nullptr;
    std::string sqlcmd;
    std::string db_name = name + ".db";
    rc = sqlite3_open(db_name.c_str(), &db);

    if (rc != SQLITE_OK)
    {
        std::cout << "SQL error: " << zErrMsg << std::endl;
        std::string m(zErrMsg);
        helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "WriteCampaignData",
                                             "SQL error on writing records:");
        sqlite3_free(zErrMsg);
    }
    sqlcmd = "CREATE TABLE if not exists bpfiles (name PRIMARY KEY);";
    rc = sqlite3_exec(db, sqlcmd.c_str(), 0, 0, &zErrMsg);
    if (rc != SQLITE_OK)
    {
        std::cout << "SQL error: " << zErrMsg << std::endl;
        std::string m(zErrMsg);
        helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "WriteCampaignData",
                                             "SQL error on writing records:");
        sqlite3_free(zErrMsg);
    }

    for (auto &r : cmap)
    {
        sqlcmd = "INSERT OR IGNORE INTO bpfiles (name)\n";
        sqlcmd += "VALUES('" + r.first + "');";
        rc = sqlite3_exec(db, sqlcmd.c_str(), 0, 0, &zErrMsg);
        if (rc != SQLITE_OK)
        {
            std::cout << "SQL error: " << zErrMsg << std::endl;
            std::string m(zErrMsg);
            helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "WriteCampaignData",
                                                 "SQL error on writing records:");
            sqlite3_free(zErrMsg);
        }
    }

    sqlite3_close(db);

    return 0;
}

CampaignManager::CampaignManager(adios2::helper::Comm &comm) { m_WriterRank = comm.Rank(); }

CampaignManager::~CampaignManager()
{
    if (m_Opened)
    {
        Close();
    }
}

void CampaignManager::Open(const std::string &name, const UserOptions &options)
{
    const UserOptions::Campaign &opts = options.campaign;
    m_Options.active = opts.active;
    m_Options.hostname = opts.hostname;
    m_Options.campaignstorepath = opts.campaignstorepath;
    m_Options.cachepath = opts.cachepath;
    m_Options.verbose = opts.verbose;

    m_Name = m_CampaignDir + PathSeparator + name + "_" + std::to_string(m_WriterRank);
    if (m_Options.verbose > 0)
    {
        std::cout << "Campaign Manager " << m_WriterRank << " Open(" << m_Name << ")\n";
    }
    m_Opened = true;
}

void CampaignManager::Record(const std::string &name, const size_t step, const double time)
{
    if (m_Options.verbose > 0)
    {
        std::cout << "Campaign Manager " << m_WriterRank << "   Record name = " << name
                  << " step = " << step << " time = " << time << "\n";
    }
    auto r = cmap.find(name);
    if (r != cmap.end())
    {
        // update record
        size_t last_step = r->second.steps.back();
        size_t delta_step = step - last_step;
        double last_time = r->second.times.back();
        double delta_time = time - last_time;
        auto nsteps = r->second.steps.size();
        if (nsteps == 1)
        {
            r->second.delta_step = delta_step;
            r->second.delta_time = delta_time;
        }
        else
        {
            size_t old_delta_step = r->second.steps.back() - r->second.steps.rbegin()[1];
            if (old_delta_step != delta_step)
            {
                r->second.delta_step = 0;
                r->second.delta_time = 0.0;
                r->second.varying_deltas = true;
            }
        }
        r->second.steps.push_back(step);
        r->second.times.push_back(time);
    }
    else
    {
        // new entry
        CampaignRecord r(step, time);
        cmap.emplace(name, r);
    }
}

void CampaignManager::Close()
{
    if (!cmap.empty())
    {
        helper::CreateDirectory(m_CampaignDir);
        CMapToSqlite(cmap, m_WriterRank, m_Name);
    }
    m_Opened = false;
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
