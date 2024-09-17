/*
 * adiosStream.h
 *
 *  Created on: Nov 2018
 *      Author: Norbert Podhorszki
 */

#ifndef ADIOSSTREAM_H
#define ADIOSSTREAM_H

#include <cstdio>

#include "adios2.h"
#include "stream.h"

class adiosStream : public Stream
{
public:
    adios2::Engine engine;
    adios2::IO io;
    bool hasIOTimer;
    size_t appID;

    adiosStream(const std::string &streamName, adios2::IO &io, const adios2::Mode mode,
                MPI_Comm comm, bool iotimer, size_t appid);
    ~adiosStream();
    void Write(CommandWrite *cmdW, Config &cfg, const Settings &settings, size_t step);
    adios2::StepStatus Read(CommandRead *cmdR, Config &cfg, const Settings &settings, size_t step);
    void Close();

private:
    MPI_Comm comm;
    FILE *perfLogFP;
    double openTime;
    void defineADIOSArray(const std::shared_ptr<VariableInfo> ov);
    void putADIOSArray(const std::shared_ptr<VariableInfo> ov);
    void getADIOSArray(std::shared_ptr<VariableInfo> ov);
    adios2::StepStatus readADIOS(CommandRead *cmdR, Config &cfg, const Settings &settings,
                                 size_t step);
    void writeADIOS(CommandWrite *cmdW, Config &cfg, const Settings &settings, size_t step);
};

#endif /* ADIOSSTREAM_H */
