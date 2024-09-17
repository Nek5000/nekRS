/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO.h
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#ifndef IO_H_
#define IO_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include <mpi.h>

class IO
{
public:
    IO(const Settings &s, MPI_Comm comm);
    ~IO();
    void write(int step, const HeatTransfer &ht, const Settings &s, MPI_Comm comm);

private:
    std::string m_outputfilename;

    // Generate a file name from the outputfile string and the arguments
    // default is add suffix if not already there
    // if rank and step is specified, it will create a unique file name for that
    // rank and step
    std::string MakeFilename(const std::string &outputfile, const std::string &suffix,
                             int rank = -1, int step = -1)
    {
        std::string name;
        const size_t ss = outputfile.size();
        if (rank == -1 && step == -1)
        {
            // add suffix if not present already
            name = outputfile;

            // if it is a short filename, add directly
            if (ss <= suffix.size())
            {
                name += suffix;
            }

            // otherwise check whether suffix is already there
            if ((ss > suffix.size()) && outputfile.find(suffix) != ss - suffix.size())
            {
                name += suffix;
            }
        }
        else
        {
            // we need a unique name here
            name = outputfile;
            if ((ss > suffix.size()) && outputfile.find(suffix) == ss - suffix.size())
            {
                name = outputfile.substr(0, ss - suffix.size());
            }
            if (rank >= 0)
            {
                std::string rs = std::to_string(rank);
                name += "." + rs;
            }
            if (step >= 0)
            {
                std::string ts = std::to_string(step);
                name += "." + ts;
            }
            name += suffix;
        }
        return name;
    }
};

#endif /* IO_H_ */
