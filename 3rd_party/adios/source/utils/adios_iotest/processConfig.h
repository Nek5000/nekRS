/*
 * processConfig.h
 *
 *  Created on: Oct 2018
 *      Author: Norbert Podhorszki
 */

#ifndef PROCESS_CONFIG_H
#define PROCESS_CONFIG_H

#include <map>
#include <string>
#include <vector>

#include "adios2.h"

#include "settings.h"

struct VariableInfo
{
    std::string name;
    std::string type;
    adios2::ShapeID shapeID;
    size_t ndim;
    adios2::Dims shape;
    adios2::Dims decomp;
    adios2::Dims start;
    adios2::Dims count;
    size_t elemsize;
    size_t datasize;
    std::vector<char> data;
    bool readFromInput;
};

enum Operation
{
    Sleep,
    Busy,
    Write,
    Read
};

class Command
{
public:
    Operation op;
    std::string conditionalStream;
    Command(Operation operation);
    virtual ~Command() = 0;
};

class CommandSleep : public Command
{
public:
    const size_t sleepTime_us = 0; // in microseconds
    CommandSleep(size_t time);
    ~CommandSleep();
};

class CommandBusy : public Command
{
public:
    const size_t cycles = 1;
    const size_t busyTime_us = 0; // in microseconds per cycle
    CommandBusy(size_t cycles, size_t time);
    ~CommandBusy();
};

class CommandWrite : public Command
{
public:
    const std::string streamName;
    const std::string groupName;
    std::vector<std::shared_ptr<VariableInfo>> variables;
    CommandWrite(std::string stream, std::string group);
    ~CommandWrite();
};

class CommandRead : public Command
{
public:
    const adios2::StepMode stepMode;
    const std::string streamName;
    const std::string groupName;
    const float timeout_sec;
    std::vector<std::shared_ptr<VariableInfo>> variables;
    CommandRead(std::string stream, std::string group, const float timeoutSec = -1.0);
    ~CommandRead();
};

struct Config
{
    size_t nSteps = 1;
    // list of input streams that we loop over instead of nSteps
    std::map<std::string, bool> stepOverStreams;
    // groupName, list of variables to preserve user defined order
    std::map<std::string, std::vector<std::shared_ptr<VariableInfo>>> groupVariableListMap;
    // same group/variables but in an ordered map for finding
    // a particular variable
    std::map<std::string, std::map<std::string, std::shared_ptr<VariableInfo>>> groupVariablesMap;
    // appID, list of commands
    std::vector<std::shared_ptr<Command>> commands;
    // Read streams status flag for supporting conditionals
    std::map<std::string, adios2::StepStatus> condMap; // stream name
};

const std::vector<std::pair<std::string, size_t>> supportedTypes = {
    {"double", sizeof(double)}, {"float", sizeof(float)}, {"int", sizeof(int)}};

Config processConfig(const Settings &settings, size_t *currentConfigLineNumber);

#endif /* PROCESS_CONFIG_H */
