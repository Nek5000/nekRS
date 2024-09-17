/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Read local arrays from multiple processors using BlocksInfo and
 * block-selection.
 *
 * If one cannot or does not want to organize arrays present on each process
 * as one global array, still one can write them out with the same name.
 * Reading, however, needs to be handled differently: each process' array has
 * to be read separately, using Block selections. The size of each process
 * block should be discovered by the reading application by inquiring per-block
 * size information of the variable, and allocate memory for reading
 * accordingly.
 *
 * In this example we read v0, v1, v2 and v3, in 5 output steps
 * from localArray_write, where
 * v0 has the same size on every process at every step
 * v1 has different size on each process but fixed over time
 * v2 has different size on each process and that is changing over time
 * v3 is like v2 but also the number of processes writing it changes over time
 *
 * The reading method shown here works for global arrays just as well, since
 * global arrays are nothing but local arrays with extra metadata to present
 * them as global arrays. The only difference is that the block.Start array
 * tells the offset of the block in the global space, while it's empty in
 * this example.
 *
 * Created on: Aug 20, 2019
 *      Author: Norbert Podhorszki <pnorbert@ornl.gov>
 */

#include <algorithm> //std::for_each
#include <array>
#include <chrono>
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <stdexcept> //std::invalid_argument std::exception
#include <string>
#include <thread>
#include <vector>

#include <adios2.h>

std::string DimsToString(const adios2::Dims &dims)
{
    std::string s = "\"";
    for (size_t i = 0; i < dims.size(); i++)
    {
        if (i > 0)
        {
            s += ", ";
        }
        s += std::to_string(dims[i]);
    }
    s += "\"";
    return s;
}

void ReadVariable(const std::string &name, adios2::IO &io, adios2::Engine &reader, size_t step)
{
    adios2::Variable<double> variable = io.InquireVariable<double>(name);

    if (variable)
    {
        auto blocksInfo = reader.BlocksInfo(variable, step);

        std::cout << "    " << name << " has " << blocksInfo.size() << " blocks in step " << step
                  << std::endl;

        // create a data vector for each block
        std::vector<std::vector<double>> dataSet;
        dataSet.resize(blocksInfo.size());

        // schedule a read operation for each block separately
        int i = 0;
        for (auto &info : blocksInfo)
        {
            variable.SetBlockSelection(info.BlockID);
            reader.Get<double>(variable, dataSet[i], adios2::Mode::Deferred);
            ++i;
        }

        // Read in all blocks at once now
        reader.PerformGets();
        // data vectors now are filled with data

        i = 0;
        for (const auto &info : blocksInfo)
        {
            std::cout << "        block " << info.BlockID << " size = " << DimsToString(info.Count)
                      << " offset = " << DimsToString(info.Start) << " : ";

            for (const auto datum : dataSet[i])
            {
                std::cout << datum << " ";
            }
            std::cout << std::endl;
            ++i;
        }
    }
    else
    {
        std::cout << "    Variable " << name << " not found in step " << step << std::endl;
    }
}

int main(int argc, char *argv[])
{
    try
    {
        adios2::ADIOS adios;

        /*** IO class object: settings and factory of Settings: Variables,
         * Parameters, Transports, and Execution: Engines
         * Inline uses single IO for write/read */
        adios2::IO io = adios.DeclareIO("Input");

        io.SetParameters({{"verbose", "4"}});

        adios2::Engine reader = io.Open("localArray.bp", adios2::Mode::Read);

        while (true)
        {

            // Begin step
            adios2::StepStatus read_status = reader.BeginStep(adios2::StepMode::Read, 10.0f);
            if (read_status == adios2::StepStatus::NotReady)
            {
                // std::cout << "Stream not ready yet. Waiting...\n";
                std::this_thread::sleep_for(std::chrono::milliseconds(1000));
                continue;
            }
            else if (read_status != adios2::StepStatus::OK)
            {
                break;
            }

            size_t step = reader.CurrentStep();
            std::cout << "Process step " << step << ": " << std::endl;

            ReadVariable("v0", io, reader, step);
            ReadVariable("v1", io, reader, step);
            ReadVariable("v2", io, reader, step);
            ReadVariable("v3", io, reader, step);

            reader.EndStep();
        }

        reader.Close();
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM\n";
        std::cout << e.what() << "\n";
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM\n";
        std::cout << e.what() << "\n";
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank\n";
        std::cout << e.what() << "\n";
    }
}
