/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * inlineMixedLang.cpp example borrowed from bpFWriteCRead example, but using
 * the inline engine. inlineMixedLang.f90 creates data and uses inline writer
 * for making that data available in C++ functions that reads the data.
 *
 *  Created on: Oct 20, 2021
 *      Author: Caitlin Ross caitlin.ross@kitware.com
 */

#include <iostream>

#include "InlineExampleFC.h"
#include <adios2.h>

adios2::IO *adiosIO;
adios2::Engine reader;

extern "C" {

void FC_GLOBAL(open_reader, OPEN_READER)(adios2::IO *io)
{
    adiosIO = io;
    reader = adiosIO->Open("reader", adios2::Mode::Read);
}

void FC_GLOBAL(analyze_data, ANALYZE_DATA)()
{
    // begin step on the reader
    reader.BeginStep();

    // grab the desired variable
    auto u = adiosIO->InquireVariable<float>("data2D");
    if (!u)
    {
        std::cerr << "variable data2D not found!" << std::endl;
        return;
    }
    auto shape = u.Shape();

    // get the data pointer and do some stuff with it
    float *data = nullptr;
    reader.Get(u, &data);
    std::cout << "c++ output: " << std::endl;
    for (size_t i = 0; i < shape[1]; i++)
    {
        for (size_t j = 0; j < shape[0]; j++)
        {
            std::cout << data[shape[0] * i + j] << "\t";
        }
        std::cout << std::endl;
    }

    // end step for the reader, signaling that we're done using the data
    reader.EndStep();
}

void FC_GLOBAL(close_reader, CLOSE_READER)() { reader.Close(); }
}
