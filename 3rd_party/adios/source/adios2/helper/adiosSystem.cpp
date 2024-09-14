/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosSystem.cpp implementation of adiosSystem.h functions
 *
 *  Created on: May 17, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include "adiosSystem.h"
#include <chrono> //system_clock, now
#include <ctime>
#include <stdexcept> // std::runtime_error, std::exception
#include <system_error>
#include <thread>

#ifndef _WIN32
#include <sys/resource.h> // getrlimits, setrlimits
#include <sys/time.h>
#endif

#include <adios2sys/SystemTools.hxx>

#include "adios2/common/ADIOSTypes.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosString.h"

// needed by IsHDF5File()
#include "adios2/core/IO.h"
#include "adios2/toolkit/transportman/TransportMan.h"
#include <cstring>

// remove ctime warning on Windows
#ifdef _WIN32
#pragma warning(disable : 4996) // ctime warning
#endif

namespace adios2
{
namespace helper
{

bool CreateDirectory(const std::string &fullPath) noexcept
{
    return static_cast<bool>(adios2sys::SystemTools::MakeDirectory(fullPath));
}

bool IsLittleEndian() noexcept
{
    uint16_t hexa = 0x1234;
    return *reinterpret_cast<uint8_t *>(&hexa) != 0x12; // NOLINT
}

std::string LocalTimeDate() noexcept
{
    struct tm now_tm;
    char buf[30];

    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

#ifdef _WIN32
    localtime_s(&now_tm, &now);
#else
    localtime_r(&now, &now_tm);
#endif
    strftime(buf, sizeof(buf), "%a %b %d %H:%M:%S %Y\n", &now_tm);

    return std::string(buf);
}

bool IsRowMajor(const std::string hostLanguage) noexcept
{
    bool isRowMajor = true;

    if (hostLanguage == "Fortran" || hostLanguage == "R" || hostLanguage == "Matlab")
    {
        isRowMajor = false;
    }

    return isRowMajor;
}

bool IsZeroIndexed(const std::string hostLanguage) noexcept
{
    bool isZeroIndexed = true;

    if (hostLanguage == "Fortran" || hostLanguage == "R")
    {
        isZeroIndexed = false;
    }
    return isZeroIndexed;
}

int ExceptionToError(const std::string &function)
{
    try
    {
        throw;
    }
    catch (std::invalid_argument &e)
    {
        helper::Log("Helper", "adiosSystem", "ExceptionToError", function + ": " + e.what(),
                    helper::FATALERROR);
        return 1;
    }
    catch (std::system_error &e)
    {
        helper::Log("Helper", "adiosSystem", "ExceptionToError", function + ": " + e.what(),
                    helper::FATALERROR);
        return 2;
    }
    catch (std::runtime_error &e)
    {
        helper::Log("Helper", "adiosSystem", "ExceptionToError", function + ": " + e.what(),
                    helper::FATALERROR);
        return 3;
    }
    catch (std::exception &e)
    {
        helper::Log("Helper", "adiosSystem", "ExceptionToError", function + ": " + e.what(),
                    helper::FATALERROR);
        return 4;
    }
}

bool IsHDF5File(const std::string &name, core::IO &io, helper::Comm &comm,
                const std::vector<Params> &transportsParameters) noexcept
{
    bool isHDF5 = false;
    if (!comm.Rank())
    {
        try
        {
            transportman::TransportMan tm(io, comm);
            if (transportsParameters.empty())
            {
                std::vector<Params> defaultTransportParameters(1);
                defaultTransportParameters[0]["transport"] = "File";
                tm.OpenFiles({name}, adios2::Mode::Read, defaultTransportParameters, false);
            }
            else
            {
                tm.OpenFiles({name}, adios2::Mode::Read, transportsParameters, false);
            }
            const unsigned char HDF5Header[8] = {137, 72, 68, 70, 13, 10, 26, 10};
            if (tm.GetFileSize(0) >= 8)
            {
                char header[8];
                tm.ReadFile(header, 8, 0);
                tm.CloseFiles();
                isHDF5 = !std::memcmp(header, HDF5Header, 8);
            }
        }
        catch (std::ios_base::failure &)
        {
            isHDF5 = false;
        }
    }
    size_t flag = (isHDF5 ? 1 : 0);
    flag = comm.BroadcastValue(flag);
    return (flag == 1);
}

char BPVersion(const std::string &name, helper::Comm &comm,
               const std::vector<Params> &transportsParameters) noexcept
{
    char version = '4'; // default result
    if (!comm.Rank())
    {
        std::string mmdFileName = name + PathSeparator + "mmd.0";
        if (adios2sys::SystemTools::PathExists(mmdFileName))
        {
            version = '5';
        }
        else
        {
            version = '4';
        }
    }
    version = comm.BroadcastValue(version);
    return version;
}

unsigned int NumHardwareThreadsPerNode() { return std::thread::hardware_concurrency(); }

size_t RaiseLimitNoFile()
{
#ifdef _WIN32
    return _setmaxstdio(8192);
#else
    static size_t raisedLimit = 0;
    static bool firstCallRaiseLimit = true;

    if (firstCallRaiseLimit)
    {
        struct rlimit limit;
        errno = 0;
        int err = getrlimit(RLIMIT_NOFILE, &limit);
        raisedLimit = limit.rlim_cur;
        if (!err)
        {
            /*std::cout
                << "adios2::helper::RaiseLimitNoFile() found limits soft = "
                << limit.rlim_cur << " hard = " << limit.rlim_max <<
               std::endl;*/
            if (limit.rlim_cur < limit.rlim_max)
            {
                limit.rlim_cur = limit.rlim_max;
                err = setrlimit(RLIMIT_NOFILE, &limit);
                if (!err)
                {
                    getrlimit(RLIMIT_NOFILE, &limit);
                    raisedLimit = limit.rlim_cur;
                    /*std::cout << "adios2::helper::RaiseLimitNoFile() set "
                                 "limits soft = "
                              << limit.rlim_cur << " hard = " << limit.rlim_max
                              << std::endl;*/
                }
            }
        }

        if (err)
        {
            std::cerr << "adios2::helper::RaiseLimitNoFile(soft=" << limit.rlim_cur
                      << ", hard=" << limit.rlim_max << ") failed with error code " << errno << ": "
                      << strerror(errno) << std::endl;
        }
        firstCallRaiseLimit = false;
    }
    return raisedLimit;
#endif
}

} // end namespace helper
} // end namespace adios2
