/*
 * stream.h
 *
 *  Created on: Nov 2018
 *      Author: Norbert Podhorszki
 */

#ifndef STREAM_H
#define STREAM_H

#include "adios2.h"
#include "ioGroup.h"
#include "processConfig.h"
#include "settings.h"

#include <string>

class Stream
{
public:
    const std::string streamName;
    adios2::Mode mode;
    Stream(const std::string &streamName, const adios2::Mode mode);
    virtual ~Stream() = 0;
    virtual void Write(CommandWrite *cmdW, Config &cfg, const Settings &settings, size_t step) = 0;
    virtual adios2::StepStatus Read(CommandRead *cmdR, Config &cfg, const Settings &settings,
                                    size_t step) = 0;
    virtual void Close() = 0;

protected:
    void fillArray(std::shared_ptr<VariableInfo> ov, double value);
};

std::shared_ptr<Stream> openStream(const std::string &streamName, std::shared_ptr<ioGroup> iogroup,
                                   const adios2::Mode mode, IOLib iolib, MPI_Comm comm,
                                   bool iotimer, size_t appid);

#endif /* STREAM_H */
