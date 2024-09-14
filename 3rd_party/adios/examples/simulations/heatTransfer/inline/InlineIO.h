/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * InlineIO.h
 *
 *  Created on: May 2020
 *      Author: Caitlin Ross
 */

#ifndef InlineIO_H_
#define InlineIO_H_

#include "../write/HeatTransfer.h"
#include "../write/Settings.h"

#include <mpi.h>

#include "adios2.h"

class InlineIO
{
public:
    InlineIO(const Settings &s, MPI_Comm comm);
    ~InlineIO();
    void write(const HeatTransfer &ht);
    const double *read(bool firstStep);

private:
    adios2::ADIOS ad;
    adios2::IO inlineIO;
    adios2::Engine inlineWriter;
    adios2::Engine inlineReader;

    adios2::Variable<double> varT;

    // need this so data on the write side doesn't go out
    // of scope before the reader can use it
    std::vector<double> v;
};

#endif /* InlineIO_H_ */
