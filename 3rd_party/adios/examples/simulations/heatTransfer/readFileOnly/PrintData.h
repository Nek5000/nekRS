/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * PrintData.h
 *
 *  Created on: Apr 2017
 *      Author: Norbert Podhorszki
 */

#ifndef PRINTDATA_H_
#define PRINTDATA_H_

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

template <class T>
void printData(double *xy, T *size, T *offset, int rank, size_t steps)
{
    std::ofstream myfile;
    std::string filename = "data." + std::to_string(rank);
    myfile.open(filename);
    double *data = xy;
    uint64_t nelems = size[0] * size[1];
    for (size_t step = 0; step < steps; step++)
    {
        myfile << "rank=" << rank << " size=" << size[0] << "x" << size[1]
               << " offsets=" << offset[0] << ":" << offset[1] << " step=" << step << std::endl;

        myfile << " time   row   columns " << offset[1] << "..." << offset[1] + size[1] - 1
               << std::endl;
        myfile << "        ";
        for (int j = 0; j < static_cast<int>(size[1]); j++)
        {
            myfile << std::setw(9) << offset[1] + j;
        }
        myfile << std::endl;
        myfile << "------------------------------------------------------------"
                  "--\n";
        for (int i = 0; i < static_cast<int>(size[0]); i++)
        {
            myfile << std::setw(5) << step << std::setw(5) << offset[0] + i;
            for (int j = 0; j < static_cast<int>(size[1]); j++)
            {
                myfile << std::setw(9) << std::setprecision(2) << data[i * size[1] + j];
            }
            myfile << std::endl;
        }
        data += nelems;
    }
    myfile.close();
}

#endif /* PRINTDATA_H_ */
