/*
 * processConfig.cpp
 *
 *  Created on: Oct 2018
 *      Author: Norbert Podhorszki
 */

#include <algorithm>
#include <chrono>
#include <errno.h>
#include <float.h> // FLT_MAX
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "decomp.h"
#include "processConfig.h"

#include "adios2/helper/adiosLog.h"

Command::Command(Operation operation) : op(operation) {}
Command::~Command() {}

CommandSleep::CommandSleep(size_t time) : Command(Operation::Sleep), sleepTime_us(time) {}
CommandSleep::~CommandSleep() {}

CommandBusy::CommandBusy(size_t cycles, size_t time)
: Command(Operation::Busy), cycles(cycles), busyTime_us(time)
{
}
CommandBusy::~CommandBusy() {}

CommandWrite::CommandWrite(std::string stream, std::string group)
: Command(Operation::Write), streamName(stream), groupName(group)
{
}
CommandWrite::~CommandWrite() {}

CommandRead::CommandRead(std::string stream, std::string group, const float timeoutSec)
: Command(Operation::Read), stepMode(adios2::StepMode::Read), streamName(stream), groupName(group),
  timeout_sec(timeoutSec)
{
}
CommandRead::~CommandRead() {}

std::vector<std::string> FileToLines(std::ifstream &configfile)
{
    std::vector<std::string> lines;
    for (std::string line; std::getline(configfile, line); /**/)
    {
        lines.push_back(line);
    }
    return lines;
}

std::vector<std::string> LineToWords(const std::string &line)
{
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
              back_inserter(tokens));
    return tokens;
}

bool isComment(std::string &s)
{
    bool comment = false;
    if (!s.compare(0, 1, "#") || !s.compare(0, 1, "%"))
    {
        comment = true;
    }
    return comment;
}

size_t stringToSizet(std::vector<std::string> &words, size_t pos, std::string lineID)
{
    if (words.size() < pos + 1)
    {
        adios2::helper::Throw<std::invalid_argument>(
            "Utils::adios_iotest", "processConfig", "stringToSizet",
            "Line for " + lineID + " is invalid. Missing value at word position " +
                std::to_string(pos + 1));
    }

    char *end;
    errno = 0;
    size_t n = static_cast<size_t>(std::strtoull(words[pos].c_str(), &end, 10));
    if (end[0] || errno == ERANGE)
    {
        adios2::helper::Throw<std::invalid_argument>(
            "Utils::adios_iotest", "processConfig", "stringToSizet",
            "Invalid value given for " + lineID + ": " + words[pos]);
    }
    return n;
}

double stringToDouble(std::vector<std::string> &words, size_t pos, std::string lineID)
{
    if (words.size() < pos + 1)
    {
        adios2::helper::Throw<std::invalid_argument>(
            "Utils::adios_iotest", "processConfig", "stringToDouble",
            "Line for " + lineID + " is invalid. Missing floating point value at word position " +
                std::to_string(pos + 1));
    }

    char *end;
    errno = 0;
    double d = static_cast<double>(std::strtod(words[pos].c_str(), &end));
    if (end[0] || errno == ERANGE)
    {
        adios2::helper::Throw<std::invalid_argument>(
            "Utils::adios_iotest", "processConfig", "stringToDouble",
            "Invalid floating point value given for " + lineID + ": " + words[pos]);
    }
    return d;
}

void processSteps(std::vector<std::string> &words, Config &cfg)
{
    if (words.size() >= 3)
    {
        std::string w(words[1]);
        std::transform(w.begin(), w.end(), w.begin(), ::tolower);
        if (w == "over" && !isComment(words[2]))
        {
            cfg.nSteps = 0;
            for (size_t i = 2; words.size() > i; ++i)
            {
                if (isComment(words[i]))
                {
                    break;
                }
                cfg.stepOverStreams[words[i]] = true;
            }
            return;
        }
    }
    cfg.nSteps = stringToSizet(words, 1, "steps");
}

void PrintDims(const adios2::Dims &dims) noexcept
{
    std::cout << "{";
    for (size_t i = 0; i < dims.size(); i++)
    {
        std::cout << dims[i];
        if (i < dims.size() - 1)
            std::cout << ",";
    }
    std::cout << "}";
}

std::string DimsToString(const adios2::Dims &dims) noexcept
{
    std::string s = "{";
    for (size_t i = 0; i < dims.size(); i++)
    {
        s += std::to_string(dims[i]);
        if (i < dims.size() - 1)
            s += ",";
    }
    s += "}";
    return s;
}

size_t processDecomp(std::string &word, const Settings &settings, std::string decompID)
{
    size_t decomp = 1;
    std::string w(word);
    std::transform(w.begin(), w.end(), w.begin(), ::toupper);
    for (size_t i = 0; i < word.size(); i++)
    {
        char c = w[i];
        if (c == 'X')
        {
            decomp *= settings.processDecomp[0];
        }
        else if (c == 'Y')
        {
            decomp *= settings.processDecomp[1];
        }
        else if (c == 'Z')
        {
            decomp *= settings.processDecomp[2];
        }
        else if (c == 'V')
        {
            decomp *= settings.processDecomp[3];
        }
        else if (c == 'W')
        {
            decomp *= settings.processDecomp[4];
        }
        else if (c == '1')
        {
            decomp *= 1;
        }
        else
        {
            adios2::helper::Throw<std::invalid_argument>(
                "Utils::adios_iotest", "processConfig", "processDecomp",
                "Invalid identifier '" + std::string(1, c) + "' for " + decompID +
                    " in character position " + std::to_string(i + 1) +
                    ". Only accepted characters are XYZVW and 1");
        }
    }
    return decomp;
}

size_t getTypeSize(std::string &type)
{
    for (const auto &t : supportedTypes)
    {
        if (t.first == type)
        {
            return t.second;
        }
    }
    adios2::helper::Throw<std::invalid_argument>("Utils::adios_iotest", "processConfig",
                                                 "getTypeSize", "Type '" + type + "' is invalid. ");
    return 0;
}

VariableInfo processArray(std::vector<std::string> &words, const Settings &settings)
{
    if (words.size() < 4)
    {
        adios2::helper::Throw<std::invalid_argument>(
            "Utils::adios_iotest", "processConfig", "processArray",
            "Line for array definition is invalid. There must be at least 4 "
            "words in the line (array type name ndim)");
    }
    VariableInfo ov;
    ov.shapeID = adios2::ShapeID::GlobalArray;
    ov.type = words[1];
    ov.elemsize = getTypeSize(ov.type);
    ov.name = words[2];
    ov.ndim = stringToSizet(words, 3, "number of dimensions of array " + ov.name);
    ov.readFromInput = false;

    if (words.size() < 4 + 2 * ov.ndim)
    {
        adios2::helper::Throw<std::invalid_argument>(
            "Utils::adios_iotest", "processConfig", "processArray",
            "Line for array definition is invalid. There must be at least 4 + "
            "2*N words where N is the 4th word ndim in the line (array type "
            "name ndim dim1 ... dimN decomp1 ... decompN)");
    }

    for (size_t i = 0; i < ov.ndim; i++)
    {
        if (settings.isStrongScaling)
        {
            ov.shape.push_back(stringToSizet(words, 4 + i, "dimension " + std::to_string(i + 1)));
        }
        else
        {
            ov.count.push_back(stringToSizet(words, 4 + i, "dimension " + std::to_string(i + 1)));
        }
    }

    size_t nprocDecomp = 1;
    for (size_t i = 0; i < ov.ndim; i++)
    {
        size_t d = processDecomp(words[4 + ov.ndim + i], settings,
                                 "decomposition " + std::to_string(i + 1));
        ov.decomp.push_back(d);
        nprocDecomp *= d;
    }
    if (nprocDecomp != settings.nProc)
    {
        adios2::helper::Throw<std::invalid_argument>(
            "Utils::adios_iotest", "processConfig", "processArray",
            "Invalid decomposition for array '" + ov.name +
                "'. The product of the decompositions (here " + std::to_string(nprocDecomp) +
                ") must equal the number of processes (here " + std::to_string(settings.nProc) +
                ")");
    }
    return ov;
}

void printConfig(const Config &cfg)
{
    std::cout << "\nConfig: \n"
              << "    nSteps    =  " << cfg.nSteps
              << "    nGroups   = " << cfg.groupVariableListMap.size()
              << "    nCommands = " << cfg.commands.size() << std::endl;
    for (const auto &mapIt : cfg.groupVariableListMap)
    {
        std::cout << "    Group " << mapIt.first << ":" << std::endl;
        for (const auto &vi : mapIt.second)
        {
            std::cout << "        " << vi->type << "  " << vi->name << DimsToString(vi->shape)
                      << "  decomposed as " << DimsToString(vi->decomp) << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "    Commands :" << std::endl;
    for (const auto &cmd : cfg.commands)
    {
        if (!cmd->conditionalStream.empty())
        {
            std::cout << "        if reading " << cmd->conditionalStream
                      << " is successful then... ";
        }
        switch (cmd->op)
        {
        case Operation::Sleep: {
            auto cmdS = dynamic_cast<const CommandSleep *>(cmd.get());
            std::cout << "        Sleep for " << cmdS->sleepTime_us << " microseconds "
                      << std::endl;
            break;
        }
        case Operation::Busy: {
            auto cmdS = dynamic_cast<const CommandBusy *>(cmd.get());
            std::cout << "          Be busy for " << cmdS->cycles << " compute cycles with "
                      << cmdS->busyTime_us << " microseconds of computation each " << std::endl;
            break;
        }
        case Operation::Write: {
            auto cmdW = dynamic_cast<CommandWrite *>(cmd.get());
            std::cout << "        Write to output " << cmdW->streamName << " the group "
                      << cmdW->groupName;
            auto grpIt = cfg.groupVariablesMap.find(cmdW->groupName);
            if (cmdW->variables.size() < grpIt->second.size())
            {
                std::cout << " with selected variables:  ";
                for (const auto &v : cmdW->variables)
                {
                    std::cout << v->name << " ";
                }
            }
            std::cout << std::endl;
            break;
        }
        case Operation::Read: {
            auto cmdR = dynamic_cast<CommandRead *>(cmd.get());
            std::cout << "        Read ";
            if (cmdR->stepMode == adios2::StepMode::Read)
            {
                std::cout << "next available step from ";
            }
            else
            {
                std::cout << "latest step from ";
            }

            std::cout << cmdR->streamName << " using the group " << cmdR->groupName;
            if (!cmdR->variables.empty())
            {
                std::cout << " with selected variables:  ";
                for (const auto &v : cmdR->variables)
                {
                    std::cout << v->name << " ";
                }
            }
            std::cout << std::endl;
            break;
        }
        }
        std::cout << std::endl;
    }
}

void printVarMaps(Config &cfg, std::string &groupName)
{
    std::cout << "DEBUG: PrintVarMap group =  " << groupName << std::endl;
    auto grpIt = cfg.groupVariablesMap.find(groupName);
    // std::cout << "    varMap = " << static_cast<void *>(grpIt->second)
    //          << std::endl;

    for (auto &v : grpIt->second)
    {
        std::cout << "     variable name first = " << v.first
                  << " second->name = " << v.second->name << " type = " << v.second->type
                  << std::endl;
    }
}

void globalChecks(const Config &cfg, const Settings &settings)
{
    for (const auto &cmd : cfg.commands)
    {
        if (!cmd->conditionalStream.empty())
        {
            if (cfg.condMap.find(cmd->conditionalStream) == cfg.condMap.end())
            {
                adios2::helper::Throw<std::invalid_argument>(
                    "Utils::adios_iotest", "processConfig", "globalChecks",
                    "Name used in conditional is not a read stream: '" + cmd->conditionalStream +
                        "'");
            }
        }
    }
    for (const auto &it : cfg.stepOverStreams)
    {
        if (cfg.condMap.find(it.first) == cfg.condMap.end())
        {
            adios2::helper::Throw<std::invalid_argument>(
                "Utils::adios_iotest", "processConfig", "globalChecks",
                "Name used in step over command is not a read stream: '" + it.first + "' ");
        }
    }
}

Config processConfig(const Settings &settings, size_t *currentConfigLineNumber)
{
    unsigned int verbose0 = (settings.myRank ? 0 : settings.verbose); // only rank 0 prints info
    std::ifstream configFile(settings.configFileName);
    if (!configFile.is_open())
    {
        adios2::helper::Throw<std::invalid_argument>(
            "Utils::adios_iotest", "processConfig", "processConfig",
            settings.configFileName + " cannot be opened ");
    }
    if (verbose0)
    {
        std::cout << "Process config file " << settings.configFileName << std::endl;
    }

    Config cfg;
    std::string currentGroup;
    int currentAppId = -1;
    std::vector<std::shared_ptr<VariableInfo>> *currentVarList = nullptr;
    std::map<std::string, std::shared_ptr<VariableInfo>> *currentVarMap = nullptr;
    std::vector<std::string> lines = FileToLines(configFile);
    for (auto &line : lines)
    {
        std::string conditionalStream;
        ++*currentConfigLineNumber;
        if (verbose0 > 1)
        {
            std::cout << "config " << *currentConfigLineNumber << ": " << line << std::endl;
        }
        std::vector<std::string> words = LineToWords(line);
        if (!words.empty() && !isComment(words[0]))
        {
            std::string key(words[0]);
            std::transform(key.begin(), key.end(), key.begin(), ::tolower);
            if (key == "cond")
            {
                if (words.size() < 2)
                {
                    adios2::helper::Throw<std::invalid_argument>(
                        "Utils::adios_iotest", "processConfig", "processConfig",
                        "Line for 'cond' is invalid. Missing group name at "
                        "word position 2");
                }
                conditionalStream = words[1];

                if (words.size() < 3)
                {
                    adios2::helper::Throw<std::invalid_argument>(
                        "Utils::adios_iotest", "processConfig", "processConfig",
                        "Line for 'cond' is invalid. Missing command from word "
                        "position 3");
                }
                words.erase(words.begin(), words.begin() + 2);
                key = words[0];
                std::transform(key.begin(), key.end(), key.begin(), ::tolower);
            }

            if (key == "group")
            {
                if (words.size() < 2)
                {
                    adios2::helper::Throw<std::invalid_argument>(
                        "Utils::adios_iotest", "processConfig", "processConfig",
                        "Line for group is invalid. Missing name at word "
                        "position 2");
                }

                currentGroup = words[1];
                if (verbose0)
                {
                    std::cout << "--> New variable group: " << currentGroup << std::endl;
                }
                auto it1 = cfg.groupVariableListMap.emplace(
                    currentGroup, std::initializer_list<std::shared_ptr<VariableInfo>>{});
                currentVarList = &it1.first->second;
                std::map<std::string, std::shared_ptr<VariableInfo>> emptymap;
                auto it2 = cfg.groupVariablesMap.emplace(currentGroup, emptymap);
                currentVarMap = &it2.first->second;
            }
            else if (key == "app")
            {
                currentAppId = static_cast<int>(stringToSizet(words, 1, "app"));
                if (verbose0)
                {
                    std::cout << "--> Application ID is set to: " << currentAppId;
                    if (currentAppId != static_cast<int>(settings.appId))
                    {
                        std::cout << "  Ignore commands set for this ID" << std::endl;
                    }
                    else
                    {
                        std::cout << "  This is our app" << std::endl;
                    }
                }
            }
            else if (key == "steps")
            {
                if (currentAppId == static_cast<int>(settings.appId))
                {
                    processSteps(words, cfg);
                    if (verbose0)
                    {
                        if (cfg.nSteps)
                        {
                            std::cout << "--> Steps is set to: " << cfg.nSteps << std::endl;
                        }
                        else
                        {
                            std::cout << "--> Steps over streams: ";
                            for (const auto &it : cfg.stepOverStreams)
                            {
                                std::cout << it.first << "  ";
                            }
                            std::cout << std::endl;
                        }
                    }
                }
            }
            else if (key == "sleep")
            {
                if (currentAppId == static_cast<int>(settings.appId))
                {
                    double d = stringToDouble(words, 1, "sleep");
                    if (verbose0)
                    {
                        std::cout << "--> Command Sleep for: " << std::setprecision(7) << d
                                  << " seconds" << std::endl;
                    }
                    size_t t_us = static_cast<size_t>(d * 1000000);
                    auto cmd = std::make_shared<CommandSleep>(t_us);
                    cmd->conditionalStream = conditionalStream;
                    cfg.commands.push_back(cmd);
                }
            }
            else if (key == "busy")
            {
                if (currentAppId == static_cast<int>(settings.appId))
                {
                    size_t cycles = static_cast<int>(stringToSizet(words, 1, "busy"));
                    double d = stringToDouble(words, 2, "busy");
                    if (verbose0)
                    {
                        std::cout << "--> Command Busy for: " << cycles << " cycles with "
                                  << std::setprecision(7) << d << " seconds computation each"
                                  << std::endl;
                    }
                    size_t t_us = static_cast<size_t>(d * 1000000);
                    auto cmd = std::make_shared<CommandBusy>(cycles, t_us);
                    cmd->conditionalStream = conditionalStream;
                    cfg.commands.push_back(cmd);
                }
            }
            else if (key == "write")
            {
                if (currentAppId == static_cast<int>(settings.appId))
                {
                    if (words.size() < 3)
                    {
                        adios2::helper::Throw<std::invalid_argument>(
                            "Utils::adios_iotest", "processConfig", "processConfig",
                            "Line for 'write' is invalid. Need at least output "
                            "name and group name ");
                    }
                    std::string fileName(words[1]);
                    std::string groupName(words[2]);
                    auto grpIt = cfg.groupVariablesMap.find(groupName);
                    if (grpIt == cfg.groupVariablesMap.end())
                    {
                        adios2::helper::Throw<std::invalid_argument>(
                            "Utils::adios_iotest", "processConfig", "processConfig",
                            "Group '" + groupName + "' used in 'write' command is undefined. ");
                    }

                    if (verbose0)
                    {
                        std::cout << "--> Command Write output = " << fileName
                                  << "  group = " << groupName << std::endl;
                    }
                    auto cmd = std::make_shared<CommandWrite>(fileName, groupName);
                    cmd->conditionalStream = conditionalStream;
                    cfg.commands.push_back(cmd);

                    // parse the optional variable list
                    size_t widx = 3;
                    while (words.size() > widx && !isComment(words[widx]))
                    {
                        auto vIt = grpIt->second.find(words[widx]);
                        if (vIt == grpIt->second.end())
                        {
                            adios2::helper::Throw<std::invalid_argument>(
                                "Utils::adios_iotest", "processConfig", "processConfig",
                                "Group '" + groupName +
                                    "' used in 'write' command has no variable "
                                    "'" +
                                    words[widx] + "' defined.");
                        }
                        cmd->variables.push_back(vIt->second);
                        ++widx;
                    }

                    if (cmd->variables.empty())
                    {
                        // no variables, we copy here ALL variables in the
                        // group
                        // (in the user defined order; copy only a pointer)
                        auto vars = cfg.groupVariableListMap.find(groupName);
                        for (auto v : vars->second)
                        {
                            cmd->variables.push_back(v);
                        }
                    }
                }
            }
            else if (key == "read")
            {
                if (currentAppId == static_cast<int>(settings.appId))
                {
                    if (words.size() < 4)
                    {
                        adios2::helper::Throw<std::invalid_argument>(
                            "Utils::adios_iotest", "processConfig", "processConfig",
                            "Line for 'read' is invalid. Need at least 3 "
                            "arguments: mode, output name, group name ");
                    }
                    std::string mode(words[1]);
                    std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
                    std::string streamName(words[2]);
                    std::string groupName(words[3]);
                    if (verbose0)
                    {
                        printVarMaps(cfg, groupName);
                    }
                    auto grpIt = cfg.groupVariablesMap.find(groupName);
                    if (grpIt == cfg.groupVariablesMap.end())
                    {
                        adios2::helper::Throw<std::invalid_argument>(
                            "Utils::adios_iotest", "processConfig", "processConfig",
                            "Group '" + groupName + "' used in 'read' command is undefined. ");
                    }
                    if (mode != "next" && mode != "latest")
                    {
                        adios2::helper::Throw<std::invalid_argument>(
                            "Utils::adios_iotest", "processConfig", "processConfig",
                            "Mode (1st argument) for 'read' is invalid. It "
                            "must be either 'next' or 'latest'");
                    }

                    double d = -1.0;
                    if (words.size() >= 4)
                    {
                        // next word is timeout
                        d = stringToDouble(words, 4, "read timeout");
                        if (d < 0.0)
                        {
                            d = -1.0;
                        }
                    }

                    if (verbose0)
                    {
                        std::cout << "--> Command Read mode = " << mode << "  input = " << words[2]
                                  << "  group = " << groupName << " timeout = ";
                        if (d < 0.0)
                        {
                            std::cout << "forever";
                        }
                        else
                        {
                            std::cout << d << " seconds";
                        }
                        std::cout << std::endl;
                    }
                    auto cmd =
                        std::make_shared<CommandRead>(streamName, groupName, static_cast<float>(d));
                    cmd->conditionalStream = conditionalStream;
                    cfg.commands.push_back(cmd);
                    cfg.condMap[streamName] = adios2::StepStatus::OK;

                    // parse the optional variable list
                    size_t widx = 5;
                    while (words.size() > widx && !isComment(words[widx]))
                    {
                        auto vIt = grpIt->second.find(words[widx]);
                        if (vIt == grpIt->second.end())
                        {
                            adios2::helper::Throw<std::invalid_argument>(
                                "Utils::adios_iotest", "processConfig", "processConfig",
                                "Group '" + groupName +
                                    "' used in 'write' command has no variable "
                                    "'" +
                                    words[widx] + "' defined.");
                        }
                        if (verbose0)
                        {
                            std::cout << "       select variable = " << vIt->second->name
                                      << std::endl;
                            std::cout << "       DEBUG variable = " << vIt->second->name
                                      << " type = " << vIt->second->type
                                      << " varmap = " << static_cast<void *>(vIt->second.get())
                                      << std::endl;
                        }
                        cmd->variables.push_back(vIt->second);
                        ++widx;
                    }

                    if (cmd->variables.empty())
                    {
                        // no variables, we copy here ALL variables in the
                        // group
                        // (in the user defined order; copy only a pointer)
                        auto vars = cfg.groupVariableListMap.find(groupName);
                        for (auto v : vars->second)
                        {
                            cmd->variables.push_back(v);
                        }
                    }
                }
            }
            else if (key == "array")
            {
                // process config line and get global array info
                std::shared_ptr<VariableInfo> ovp =
                    std::make_shared<VariableInfo>(processArray(words, settings));
                ovp->datasize = ovp->elemsize;
                // Position of rank in N-dim space
                std::vector<size_t> pos(ovp->ndim);
                // Calculate rank's position in N-dim space
                decompRowMajor(ovp->ndim, settings.myRank, ovp->decomp.data(), pos.data());

                if (settings.isStrongScaling)
                {
                    // Calculate the local size and offsets based on the
                    // definition
                    for (size_t i = 0; i < ovp->ndim; ++i)
                    {
                        size_t count = ovp->shape[i] / ovp->decomp[i];
                        size_t offs = count * pos[i];
                        if (pos[i] == ovp->decomp[i] - 1 && pos[i] != 0)
                        {
                            // last process in dim(i) need to write all the rest
                            // of dimension
                            count = ovp->shape[i] - offs;
                        }
                        ovp->start.push_back(offs);
                        ovp->count.push_back(count);
                        ovp->datasize *= count;
                    }
                }
                else
                {
                    // Calculate the local size and offsets based on the
                    // definition
                    for (size_t i = 0; i < ovp->ndim; ++i)
                    {
                        size_t shape = ovp->count[i] * ovp->decomp[i];
                        size_t offs = ovp->count[i] * pos[i];
                        ovp->start.push_back(offs);
                        ovp->shape.push_back(shape);
                        ovp->datasize *= ovp->count[i];
                    }
                }

                // Postpone Allocating the data array until first use
                // ovp->data.resize(ovp->datasize);

                currentVarList->push_back(ovp);
                currentVarMap->emplace(ovp->name, ovp);

                /* DEBUG */
                if (verbose0)
                {
                    auto grpIt = cfg.groupVariablesMap.find(currentGroup);
                    auto vIt = grpIt->second.find(ovp->name);
                    std::cout << "       DEBUG variable = " << vIt->second->name
                              << " type = " << vIt->second->type
                              << " varmap = " << static_cast<void *>(vIt->second.get())
                              << std::endl;
                }
                if (settings.verbose > 2)
                {
                    std::cout << "--> rank = " << settings.myRank
                              << ": Variable array name = " << ovp->name << " type = " << ovp->type
                              << " elemsize = " << ovp->elemsize
                              << " local datasize = " << ovp->datasize
                              << " shape = " << DimsToString(ovp->shape)
                              << " start = " << DimsToString(ovp->start)
                              << " count = " << DimsToString(ovp->count) << std::endl;
                }
                else if (verbose0)
                {
                    std::cout << "--> Variable array name = " << ovp->name
                              << " type = " << ovp->type << " elemsize = " << ovp->elemsize
                              << std::endl;
                }
            }
            else if (key == "link")
            {
                if (words.size() < 3)
                {
                    adios2::helper::Throw<std::invalid_argument>(
                        "Utils::adios_iotest", "processConfig", "processConfig",
                        "Line for 'link' is invalid. Need at least 2 words: "
                        "group <group_name> ");
                }
                std::string groupName = words[2];
                auto grpIt = cfg.groupVariablesMap.find(groupName);
                if (grpIt == cfg.groupVariablesMap.end())
                {
                    adios2::helper::Throw<std::invalid_argument>(
                        "Utils::adios_iotest", "processConfig", "processConfig",
                        "Group '" + groupName + "' used in 'link' command is undefined. ");
                }

                if (verbose0)
                {
                    std::cout << "--> Link variables from group = " << groupName << ": ";
                }

                // parse the optional variable list
                size_t widx = 3;
                if (words.size() > widx && !isComment(words[widx]))
                {
                    while (words.size() > widx && !isComment(words[widx]))
                    {
                        auto vIt = grpIt->second.find(words[widx]);
                        if (vIt == grpIt->second.end())
                        {
                            adios2::helper::Throw<std::invalid_argument>(
                                "Utils::adios_iotest", "processConfig", "processConfig",
                                "Group '" + groupName +
                                    "' used in 'link' command has no variable "
                                    "'" +
                                    words[widx] + "' defined.");
                        }
                        if (verbose0)
                        {
                            std::cout << " " << vIt->second->name;
                        }
                        currentVarList->push_back(vIt->second);
                        currentVarMap->emplace(vIt->second->name, vIt->second);
                        ++widx;
                    }
                }
                else
                {
                    // no variable list, so we link ALL variables in the
                    // group (in the user defined order)
                    if (verbose0)
                    {
                        std::cout << "all variables";
                    }
                    auto vars = cfg.groupVariableListMap.find(groupName);
                    for (auto v : vars->second)
                    {
                        currentVarList->push_back(v);
                        currentVarMap->emplace(v->name, v);
                    }
                }
                if (verbose0)
                {
                    std::cout << std::endl;
                }
            }
            else
            {
                adios2::helper::Throw<std::invalid_argument>("Utils::adios_iotest", "processConfig",
                                                             "processConfig",
                                                             "Unrecognized keyword '" + key + "'.");
            }
        }
    }
    configFile.close();
    *currentConfigLineNumber = 0;
    globalChecks(cfg, settings);
    if (verbose0 > 2)
    {
        printConfig(cfg);
    }
    return cfg;
}
